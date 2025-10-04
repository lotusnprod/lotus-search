import functools
import io
import logging
import pickle
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Dict, Iterable as It, List, Tuple, Set

import requests
from rdkit import Chem, DataStructs
from rdkit.Chem import rdSubstructLibrary
from sqlalchemy import and_
from sqlalchemy.orm import aliased

from api.models import ReferenceObject, StructureObject, TaxonObject
from chemistry_helpers import fingerprint, standardize
from sdf_helpers import find_structures_bytes_ranges, mmap_file, read_selected_ranges
from storage.models import (
    Journals,
    References,
    Structures,
    StructuresDescriptors,
    TaxoNames,
    TaxoParents,
    TaxoRankNames,
    TaxoRanks,
    Triplets,
)
from storage.storage import Storage

logging.basicConfig()
log = logging.getLogger()
log.setLevel(logging.WARNING)
requests_log = logging.getLogger("requests.packages.urllib3")
requests_log.setLevel(logging.WARNING)
requests_log.propagate = True


class DataModel:
    """In-memory + SQLite accessor for LOTUS search functionalities.

    Responsibilities:
      - Load binary pre-computed structure similarity & substructure libraries.
      - Provide lazy / cached access to large SDF file and byte ranges.
      - Offer domain-specific query helpers combining SQLAlchemy queries
        and in-memory similarity/substructure search.

    The underlying SQLite dataset is treated as immutable for the life of the
    instance; therefore selective method-level caching is safe where used.
    """

    def __init__(self, path: Path = Path("./data")):
        self.db = self.load_all_data(path)
        # Preload SDF memory map + byte ranges for efficient SDF block extraction
        self.sdf = self.load_sdf_data(path)
        self.sdf_ranges = self.load_sdf_ranges(self.sdf)
        self.storage = Storage(path)
        self.taxa_name_db = self.preload_taxa()
        self.path = path

    # ------------------------------ Loading ---------------------------------
    @classmethod
    @functools.cache
    def load_all_data(cls, path: Path):
        """Load pickled data that embeds RDKit substructure libraries.

        The pickled payload stores the binary stream of two SubstructLibrary
        instances (with & without explicit hydrogens). They are reconstructed
        here to avoid re-computation. Cached at the class level to reuse
        across multiple DataModel instances (e.g. in tests).
        """
        with open(path / "database_chemo.pkl", "rb") as f:
            data = pickle.load(f)
        new_lib = rdSubstructLibrary.SubstructLibrary()
        new_lib_h = rdSubstructLibrary.SubstructLibrary()
        with io.BytesIO(data["structure_library"]) as i:
            new_lib.InitFromStream(i)
        with io.BytesIO(data["structure_library_h"]) as i:
            new_lib_h.InitFromStream(i)
        data["structure_library"] = new_lib
        data["structure_library_h"] = new_lib_h
        return data

    @classmethod
    @functools.cache
    def load_sdf_data(cls, path: Path):
        """Memory-map the main SDF file (heavy file, so mmap is efficient)."""
        mapped_sdf = mmap_file(path / "lotus.sdf")
        return mapped_sdf

    @classmethod
    @functools.cache
    def load_sdf_ranges(cls, sdf):  # type: ignore[no-untyped-def]
        """Pre-compute (and cache) byte ranges for each SDF record."""
        ranges = find_structures_bytes_ranges(sdf)
        return ranges

    # ----------------------------- Taxonomy ---------------------------------
    def get_taxon_object_from_dict_of_tids(
        self,
        tids: Iterable[int],
    ) -> dict[int, TaxonObject]:
        """Return TaxonObject mapping for provided taxon IDs.

        Returns an empty dict if none found.
        """
        with self.storage.session() as session:
            result = (
                session.query(
                    TaxoNames.id,
                    TaxoNames.name,
                )
                .filter(TaxoNames.id.in_(tids))
                .all()
            )
            if result:
                return {row.id: TaxonObject(name=row.name) for row in result}
            return {}

    @functools.lru_cache(maxsize=10000)
    def get_taxon_object_from_tid(self, tid: int) -> dict | None:  # legacy signature
        return self.get_taxon_object_from_dict_of_tids([tid])

    @functools.lru_cache(maxsize=20000)
    def get_taxa_with_name_matching(self, query: str, exact: bool = False) -> set[int]:
        """Return IDs of taxa whose names match the query.

        If exact is False a case-sensitive SQL LIKE %query% is used.
        """
        with self.storage.session() as session:
            matcher = TaxoNames.name == query if exact else TaxoNames.name.like(f"%{query}%")
            result = session.query(TaxoNames.id).filter(matcher)
            return {row[0] for row in result}

    @functools.lru_cache(maxsize=50000)
    def get_rank_name_from_wid(self, wid: int) -> str | None:
        with self.storage.session() as session:
            result = session.get(TaxoRankNames, wid)
            return None if result is None else result.name

    def resolve_taxon(self, query: str) -> Any:
        """Call Global Names resolver (best-effort, network errors swallowed)."""
        payload = {
            "nameStrings": [query],
            "withVernaculars": False,
            "withCapitalization": True,
            "withAllMatches": True,
        }
        log.debug(f"Querying GN resolver... {payload}")
        try:
            response = requests.post(
                "https://verifier.globalnames.org/api/v1/verifications",
                json=payload,
                headers={"Content-Type": "application/json"},
            )
        except Exception as e:  # pragma: no cover - network failure resilience
            log.error(f"Impossible to connect to GN resolver {e}.")
            return None
        log.debug(response.json())
        return response.json()

    @functools.lru_cache(maxsize=50000)
    def get_ranks_string(self, wid: int) -> str:
        with self.storage.session() as session:
            result = session.query(TaxoRanks.rank_id).filter(TaxoRanks.id == wid)
            rank_ids = [row[0] for row in result]
        if rank_ids:
            n_ranks = [self.get_rank_name_from_wid(int(it)) for it in rank_ids]
            ranks = " (" + ", ".join(set(n_ranks)) + ")"
        else:
            ranks = ""
        return ranks

    # --------------------------- Structureonomy -----------------------------
    @functools.cache
    def structures_set(self) -> set[int]:
        """Return a cached set of all known structure WIDs."""
        return set(self.db["structure_wid"])

    def get_structure_object_from_sid(self, sid: int) -> dict | None:  # legacy
        return self.get_structure_object_from_dict_of_sids([sid])

    def get_structure_descriptors_from_dict_of_sids(
        self,
        sids: Iterable[int],
    ) -> dict[str, list[Any]]:
        """Return descriptor name -> list of values for given structure IDs."""
        with self.storage.session() as session:
            result = (
                session.query(
                    StructuresDescriptors.descriptor_name,
                    StructuresDescriptors.descriptor_value,
                )
                .filter(StructuresDescriptors.structure_id.in_(sids))
                .distinct()
                .all()
            )
            result_dict: dict[str, list[Any]] = {}
            if result:
                for descriptor_name, descriptor_value in result:
                    result_dict.setdefault(descriptor_name, []).append(descriptor_value)
            return result_dict

    def get_structure_sdf_from_dict_of_sids(
        self,
        sids: Iterable[int],
    ) -> str:
        """Return concatenated SDF blocks for the provided structure IDs."""
        ranges = self.sdf_ranges
        mm = self.sdf
        return "".join(mm[start:end].decode() for sid in sids for (start, end) in (ranges[sid],))

    def get_structure_object_from_dict_of_sids(
        self,
        sids: Iterable[int],
        return_descriptors: bool = False,
    ) -> dict[int, StructureObject]:
        """Return StructureObject mapping (optionally with descriptors)."""
        with self.storage.session() as session:
            descriptors: dict[str, list[Any]] = {}
            if return_descriptors:
                descriptors = self.get_structure_descriptors_from_dict_of_sids(sids)
            result = (
                session.query(
                    Structures.id,
                    Structures.smiles,
                    Structures.smiles_no_stereo,
                    Structures.inchi,
                    Structures.inchi_no_stereo,
                    Structures.inchikey,
                    Structures.inchikey_no_stereo,
                    Structures.formula,
                )
                .filter(Structures.id.in_(sids))
                .all()
            )
            if result:
                return {
                    row.id: StructureObject(
                        smiles=row.smiles,
                        smiles_no_stereo=row.smiles_no_stereo,
                        inchi=row.inchi,
                        inchi_no_stereo=row.inchi_no_stereo,
                        inchikey=row.inchikey,
                        inchikey_no_stereo=row.inchikey_no_stereo,
                        formula=row.formula,
                        descriptors=descriptors,
                    )
                    for row in result
                }
            return {}

    def get_structure_with_descriptors(self, descriptors: dict) -> set[int]:
        """Return structure IDs matching descriptor min/max constraints.

        Input descriptors dict uses the pattern <DescriptorName>_min / _max.
        """
        with self.storage.session() as session:
            query = session.query(StructuresDescriptors.structure_id)
            min_descriptors: dict[str, Any] = {}
            max_descriptors: dict[str, Any] = {}
            for key, value in descriptors.items():
                descriptor_name = key[:-4]
                if key.endswith("_min"):
                    min_descriptors[descriptor_name] = value
                if key.endswith("_max"):
                    max_descriptors[descriptor_name] = value

            query_min = None
            for descriptor_name, min_value in min_descriptors.items():
                min_condition = (
                    StructuresDescriptors.descriptor_name == descriptor_name,
                    StructuresDescriptors.descriptor_value >= min_value,
                )
                query_min = query.filter(and_(*min_condition))

            query_max = None
            for descriptor_name, max_value in max_descriptors.items():
                max_condition = (
                    StructuresDescriptors.descriptor_name == descriptor_name,
                    StructuresDescriptors.descriptor_value <= max_value,
                )
                query_max = query.filter(and_(*max_condition))

            if query_min is not None and query_max is not None:
                result = query_min.intersect(query_max).all()
            elif query_min is not None:
                result = query_min.all()
            elif query_max is not None:
                result = query_max.all()
            else:
                result = query.all()
            return {row[0] for row in result}

    def get_structure_with_formula(self, formula: str) -> set[int]:
        with self.storage.session() as session:
            result = session.query(Structures.id).filter(Structures.formula == formula)
            return {row[0] for row in result}

    def structure_get_mol_fp_and_explicit(self, query: str):  # type: ignore[no-untyped-def]
        """Return (mol, fingerprint, explicit_h_present) for a SMILES query."""
        explicit_h = "[H]" in query
        p = Chem.SmilesParserParams()
        p.removeHs = not explicit_h
        mol = Chem.MolFromSmiles(query, p)
        if not explicit_h:
            mol = standardize(mol)
        fp = fingerprint(mol)
        return mol, fp, explicit_h

    def structure_search(self, query: str) -> list[tuple[int, float]]:
        """Similarity search (Tanimoto) returning list of (WID, score)."""
        mol, fp, explicit_h = self.structure_get_mol_fp_and_explicit(query)
        db = self.db["structure_sim_h_fps"] if explicit_h else self.db["structure_sim_fps"]
        scores = DataStructs.BulkTanimotoSimilarity(fp, db)
        return [
            (wid, score)
            for wid, score in zip(self.db["structure_wid"], scores, strict=False)
        ]

    def structure_search_substructure(
        self,
        query: str,
        chirality: bool = False,
    ) -> list[tuple[int, float]]:
        """Substructure search returning (WID, tanimoto_score) list.

        The tanimoto score is computed between the query fingerprint and the
        stored fingerprint for each match; original ordering preserved.
        """
        mol, fp, explicit_h = self.structure_get_mol_fp_and_explicit(query)
        if explicit_h:
            db = self.db["structure_library_h"]
            fp_db = self.db["structure_sim_h_fps"]
        else:
            db = self.db["structure_library"]
            fp_db = self.db["structure_sim_fps"]
        iids = db.GetMatches(
            mol,
            numThreads=-1,
            maxResults=-1,
            useQueryQueryMatches=True,
            useChirality=chirality,
        )
        new_keys = [self.db["structure_wid"][iid] for iid in iids]
        out: list[tuple[int, float]] = []
        for iid, wid in zip(iids, new_keys, strict=False):
            out.append((wid, DataStructs.TanimotoSimilarity(fp, fp_db[iid])))
        return out

    def structure_get_tsv_from_scores(self, wids: list[int], scores) -> str:  # type: ignore[no-untyped-def]
        """Return TSV string for similarity results (Wikidata URL, score, SMILES)."""
        out = "Wikidata link\tSimilarity\tSmiles\n"
        structure_objects = self.get_structure_object_from_dict_of_sids(wids)
        smiles_dict = {wid: structure_object.smiles for wid, structure_object in structure_objects.items()}
        for idx, score in enumerate(scores):
            wid = wids[idx]
            smiles = smiles_dict[wid]
            out += f"http://www.wikidata.org/entity/Q{wid}\t{score:.3f}\t{smiles}\n"
        return out

    # ----------------------------- Biblionomy -------------------------------
    def get_reference_with_id(self, rid: int) -> set[int]:
        with self.storage.session() as session:
            result = session.query(References.id).filter(References.id == rid)
            return {row[0] for row in result}

    def get_reference_object_from_dict_of_rids(
        self,
        rids: Iterable[int],
    ) -> dict[int, ReferenceObject]:
        with self.storage.session() as session:
            result = (
                session.query(
                    References.id,
                    References.doi,
                    References.title,
                    References.date,
                    References.journal,
                )
                .filter(References.id.in_(rids))
                .all()
            )
            journal_ids = {row.journal for row in result}
            journal_titles = {}
            if journal_ids:
                result_journal = session.query(
                    Journals.id,
                    Journals.title,
                ).filter(Journals.id.in_(journal_ids))
                journal_titles = {journal.id: journal.title for journal in result_journal}
            if result:
                return {
                    row.id: ReferenceObject(
                        doi=row.doi,
                        title=row.title,
                        date=row.date,
                        journal=journal_titles.get(row.journal),
                    )
                    for row in result
                }
            return {}

    def get_reference_object_from_rid(self, rid: int) -> dict | None:  # legacy
        return self.get_reference_object_from_dict_of_rids([rid])

    def get_references_with_doi(self, doi: str) -> set[int]:
        with self.storage.session() as session:
            result = session.query(References.id).filter(
                References.doi.like(f"%{doi}%"),
            )
            return {row[0] for row in result}

    def get_references_with_date(
        self,
        date_min: str = None,
        date_max: str = None,
    ) -> set[int]:
        with self.storage.session() as session:
            query = session.query(References.id)
            if date_min is not None and date_max is not None:
                query = query.filter(References.date.between(date_min, date_max))
            elif date_min is not None:
                query = query.filter(References.date >= date_min)
            elif date_max is not None:
                query = query.filter(References.date <= date_max)
            result = query.all()
            return {row[0] for row in result}

    def get_references_with_journal(self, journal_title: str) -> set[int]:
        with self.storage.session() as session:
            result = (
                session.query(References.id)
                .filter(
                    References.journal.in_(
                        session.query(Journals.id).filter(
                            Journals.title.like(f"%{journal_title}%"),
                        ),
                    ),
                )
                .all()
            )
            return {row[0] for row in result}

    def get_references_with_title(self, title: str) -> set[int]:
        with self.storage.session() as session:
            result = session.query(References.id).filter(
                References.title.like(f"%{title}%"),
            )
            return {row[0] for row in result}

    # Mixonomy
    # Todo, we probably want to still return that as a set
    def get_structures_of_taxon(self, tid: int, recursive: bool = True) -> set[int]:
        matching_structures = self.storage.get_generic_of_generic(
            Triplets.structure_id,
            Triplets.taxon_id,
            tid,
        )

        if recursive:
            with self.storage.session() as session:
                # if exist
                result = session.query(TaxoParents.child_id).filter(
                    TaxoParents.parent_id == tid,
                )

                for row in result:
                    for structure in self.get_structures_of_taxon(row[0]):
                        matching_structures.add(structure)
        return matching_structures

    def get_taxa_of_structure(self, sid: int) -> set[int]:
        return self.storage.get_generic_of_generic(
            Triplets.taxon_id,
            Triplets.structure_id,
            sid,
        )

    def get_structures_of_reference(self, rid: int) -> set[int]:
        return self.storage.get_generic_of_generic(
            Triplets.structure_id,
            Triplets.reference_id,
            rid,
        )

    def get_taxa_of_reference(self, rid: int) -> set[int]:
        return self.storage.get_generic_of_generic(
            Triplets.taxon_id,
            Triplets.reference_id,
            rid,
        )

    def get_references_of_structure(self, sid: int) -> set[int]:
        return self.storage.get_generic_of_generic(
            Triplets.reference_id,
            Triplets.structure_id,
            sid,
        )

    def get_references_of_taxon(self, tid: int) -> set[int]:
        return self.storage.get_generic_of_generic(
            Triplets.reference_id,
            Triplets.taxon_id,
            tid,
        )

    def get_references_of_structures(self, structures: set[int]) -> set[int]:
        return self.storage.get_generics_of_generics(
            Triplets.reference_id,
            Triplets.structure_id,
            structures,
        )

    def get_references_of_taxa(self, taxa: set[int]) -> set[int]:
        return self.storage.get_generics_of_generics(
            Triplets.reference_id,
            Triplets.taxon_id,
            taxa,
        )

    def get_structures_of_references(self, references: set[int]) -> set[int]:
        return self.storage.get_generics_of_generics(
            Triplets.structure_id,
            Triplets.reference_id,
            references,
        )

    def get_taxa_of_structures(self, structures: set[int]) -> set[int]:
        return self.storage.get_generics_of_generics(
            Triplets.taxon_id,
            Triplets.structure_id,
            structures,
        )

    def get_taxa_of_references(self, references: set[int]) -> set[int]:
        return self.storage.get_generics_of_generics(
            Triplets.taxon_id,
            Triplets.reference_id,
            references,
        )

    def get_triplets_for(
        self,
        reference_ids: set[int] | None,
        structure_ids: set[int] | None,
        taxon_ids: set[int] | None,
    ) -> set[tuple[int, int, int]]:
        return self.storage.get_triplets_for(reference_ids, structure_ids, taxon_ids)

    def get_taxon_by_id(self, tid: int) -> set[int]:
        with self.storage.session() as session:
            result = session.query(Triplets.taxon_id).filter(Triplets.taxon_id == tid)
            return {row[0] for row in result.distinct()}

    def get_taxon_children_by_id(self, tid: int) -> set[int]:
        with self.storage.session() as session:
            # Recursive query to fetch all children for the given taxon ID
            recursive_query = (
                session.query(TaxoParents.child_id)
                .filter(TaxoParents.parent_id == tid)
                .cte(name="recursive_query", recursive=True)
            )

            alias = aliased(TaxoParents)
            recursive_query = recursive_query.union_all(
                session.query(alias.child_id).filter(
                    alias.parent_id == recursive_query.c.child_id,
                ),
            )

            result = set(session.query(recursive_query).all())

            # This is working but not recursively to get all children
            # result = session.query(TaxoParents.id).filter(
            #     TaxoParents.parent_id == tid
            # )
            return {row[0] for row in result}

    def preload_taxa(self):
        with self.storage.session() as session:
            result = session.query(TaxoNames.id, TaxoNames.name).all()
            return {row[1]: row[0] for row in result}

    def get_dict_of_taxa_from_name(self, taxon_name: str) -> dict[str, int]:
        with self.storage.session() as session:
            matcher = TaxoNames.name.like(f"{taxon_name}%")
            result = session.query(TaxoNames.name, TaxoNames.id).filter(matcher)
        return {row[0]: row[1] for row in result}
