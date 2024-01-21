import functools
import io
import logging
import pickle
from collections.abc import Iterable
from pathlib import Path

import requests
from rdkit import Chem, DataStructs
from rdkit.Chem import rdSubstructLibrary

from chemistry_helpers import fingerprint, standardize
from storage.models import (
    References,
    Structures,
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
    def __new__(cls, path: Path = Path("./data")):
        # Each instance will be the same, it is all read-only
        if not hasattr(cls, "instance"):
            cls.instance = super(DataModel, cls).__new__(cls)
            cls.instance.db = cls.load_all_data(path)
            cls.instance.storage = Storage(path)
            cls.instance.path = path
        return cls.instance

    @classmethod
    @functools.lru_cache(maxsize=None)
    def load_all_data(cls, path: Path):
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

    ### Taxonomy
    def get_dict_of_tid_to_taxon_name(self, tid: Iterable[int]) -> dict[int, str]:
        with self.storage.session() as session:
            result = session.query(TaxoNames.id, TaxoNames.name).filter(
                TaxoNames.id.in_(tid)
            )
            return {row[0]: row[1] for row in result}

    def get_taxon_name_from_tid(self, tid: int) -> str | None:
        with self.storage.session() as session:
            result = session.get(TaxoNames, tid)
            if result is None:
                return None
            return result.name

    def get_taxa_with_name_matching(self, query: str, exact=False) -> set[int]:
        with self.storage.session() as session:
            if exact:
                matcher = TaxoNames.name == query
            else:
                matcher = TaxoNames.name.like(f"%{query}%")
            result = session.query(TaxoNames.id).filter(matcher)
            return {row[0] for row in result}

    def get_rank_name_from_wid(self, wid: int) -> str | None:
        with self.storage.session() as session:
            result = session.get(TaxoRankNames, wid)
            if result is None:
                return None
            return result.name

    def resolve_taxon(self, query: str) -> any:
        query = {
            "nameStrings": [query],
            "withVernaculars": False,
            "withCapitalization": True,
            "withAllMatches": True,
        }
        log.debug(f"Querying GN resolver... {query}")

        try:
            response = requests.post(
                "https://verifier.globalnames.org/api/v1/verifications",
                json=query,
                headers={"Content-Type": "application/json"},
            )
        except Exception as e:
            log.error(f"Impossible to connect to GN resolver {e}.")
            return None
        log.debug(response.json())
        return response.json()

    def get_ranks_string(self, wid: int) -> str:
        with self.storage.session() as session:
            result = session.query(TaxoRanks.rank_id).filter(TaxoRanks.id == wid)
            rank_ids = [row[0] for row in result]
        if len(rank_ids) > 0:
            n_ranks = [self.get_rank_name_from_wid(int(it)) for it in rank_ids]
            ranks = " (" + ", ".join(set(n_ranks)) + ")"
        else:
            ranks = ""
        return ranks

    ### Structureonomy
    @functools.lru_cache(maxsize=None)
    def structures_set(self) -> set[int]:
        # TODO use DB
        return set(self.db["structure_wid"])

    def get_structure_smiles_from_sid(self, sid: int) -> str | None:
        with self.storage.session() as session:
            out = session.get(Structures, sid)
            if out is None:
                return None
            return out.smiles

    def get_dict_of_sid_to_smiles(self, sids: Iterable[int]) -> dict[int, str]:
        with self.storage.session() as session:
            result = session.query(Structures.id, Structures.smiles).filter(
                Structures.id.in_(sids)
            )
            return {row[0]: row[1] for row in result}

    def structure_get_mol_fp_and_explicit(self, query: str) -> tuple[any, any, bool]:
        explicit_h = "[H]" in query
        p = Chem.SmilesParserParams()
        p.removeHs = not explicit_h
        mol = Chem.MolFromSmiles(query, p)

        if not explicit_h:
            mol = standardize(mol)

        fp = fingerprint(mol)
        return mol, fp, explicit_h

    ## COMMENT (AR): Should we rename this to structure_search_from_smiles
    ## and have same for InChI and co and then wrap them to a `structure_search`
    ## with inchi = "InChI=1S/" in query ...
    def structure_search(self, query: str) -> list[tuple[int, float]]:
        mol, fp, explicit_h = self.structure_get_mol_fp_and_explicit(query)

        if explicit_h:
            db = self.db["structure_sim_h_fps"]
        else:
            db = self.db["structure_sim_fps"]
        scores = DataStructs.BulkTanimotoSimilarity(fp, db)
        return [(wid, score) for wid, score in zip(self.db["structure_wid"], scores)]

    def structure_search_substructure(
        self, query: str, chirality: bool = False
    ) -> list[tuple[int, float]]:
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
        # TODO letting WID for now (also in data_structures) but to keep in mind
        new_keys = [self.db["structure_wid"][iid] for iid in iids]
        out = []
        for iid, wid in zip(iids, new_keys):
            out.append((wid, DataStructs.TanimotoSimilarity(fp, fp_db[iid])))
        return out

    def structure_get_tsv_from_scores(self, wids: list[int], scores) -> str:
        out = "Wikidata link\tSimilarity\tSmiles\n"
        smiles_dict = self.get_dict_of_sid_to_smiles(wids)
        for idx, score in enumerate(scores):
            wid = wids[idx]
            smiles = smiles_dict[wid]
            out += f"http://www.wikidata.org/entity/Q{wid}\t{score:.3f}\t{smiles}\n"
        return out

    ### Biblionomy
    def get_reference_with_id(self, rid: int) -> set[int]:
        with self.storage.session() as session:
            result = session.query(References.id).filter(References.id == rid)
            return {row[0] for row in result}

    def get_dict_of_rid_to_reference_doi(self, rid: Iterable[int]) -> dict[int, str]:
        with self.storage.session() as session:
            result = session.query(References.id, References.doi).filter(
                References.id.in_(rid)
            )
            return {row[0]: row[1] for row in result}

    def get_reference_doi_from_rid(self, rid: int) -> str | None:
        with self.storage.session() as session:
            result = session.get(References, rid)
            if result is None:
                return None
            return result.doi

    def get_references_with_doi(self, doi: str) -> set[int]:
        with self.storage.session() as session:
            result = session.query(References.id).filter(
                References.doi.like(f"%{doi}%")
            )
            return {row[0] for row in result}

    ### Mixonomy
    # Todo, we probably want to still return that as a set
    def get_structures_of_taxon(self, tid: int, recursive: bool = True) -> set[int]:
        matching_structures = self.storage.get_generic_of_generic(
            Triplets.structure_id, Triplets.taxon_id, tid
        )

        if recursive:
            with self.storage.session() as session:
                # if exist
                result = session.query(TaxoParents.id).filter(
                    TaxoParents.parent_id == tid
                )

                for row in result:
                    for structure in self.get_structures_of_taxon(row[0]):
                        matching_structures.add(structure)
        return matching_structures

    def get_taxa_of_structure(self, sid: int) -> set[int]:
        return self.storage.get_generic_of_generic(
            Triplets.taxon_id, Triplets.structure_id, sid
        )

    def get_structures_of_reference(self, rid: int) -> set[int]:
        return self.storage.get_generic_of_generic(
            Triplets.structure_id, Triplets.reference_id, rid
        )

    def get_taxa_of_reference(self, rid: int) -> set[int]:
        return self.storage.get_generic_of_generic(
            Triplets.taxon_id, Triplets.reference_id, rid
        )

    def get_references_of_structure(self, sid: int) -> set[int]:
        return self.storage.get_generic_of_generic(
            Triplets.reference_id, Triplets.structure_id, sid
        )

    def get_references_of_taxon(self, tid: int) -> set[int]:
        return self.storage.get_generic_of_generic(
            Triplets.reference_id, Triplets.taxon_id, tid
        )

    def get_references_of_structures(self, structures: set[int]) -> set[int]:
        return self.storage.get_generics_of_generics(
            Triplets.reference_id, Triplets.structure_id, structures
        )

    def get_references_of_taxa(self, taxa: set[int]) -> set[int]:
        return self.storage.get_generics_of_generics(
            Triplets.reference_id, Triplets.taxon_id, taxa
        )

    def get_structures_of_references(self, references: set[int]) -> set[int]:
        return self.storage.get_generics_of_generics(
            Triplets.structure_id, Triplets.reference_id, references
        )

    def get_taxa_of_structures(self, structures: set[int]) -> set[int]:
        return self.storage.get_generics_of_generics(
            Triplets.taxon_id, Triplets.structure_id, structures
        )

    def get_taxa_of_references(self, references: set[int]) -> set[int]:
        return self.storage.get_generics_of_generics(
            Triplets.taxon_id, Triplets.reference_id, references
        )

    def get_triplets_for(
        self,
        reference_ids: set[int] | None,
        structure_ids: set[int] | None,
        taxon_ids: set[int] | None,
    ) -> set[tuple[int, int, int]]:
        return self.storage.get_triplets_for(reference_ids, structure_ids, taxon_ids)

    def get_taxon_by_id(self, taxon_wid: int) -> set[int]:
        with self.storage.session() as session:
            result = session.query(Triplets.taxon_id).filter(
                Triplets.taxon_id == taxon_wid
            )
            return {row[0] for row in result.distinct()}
