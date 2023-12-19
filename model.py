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
        return cls.instance

    @classmethod
    @functools.lru_cache(maxsize=None)
    def load_all_data(cls, path: Path):
        with open(path / "database.pkl", "rb") as f:
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

    def num_taxa(self):
        return len(self.db["taxonomy_names"])

    def num_structures(self):
        return len(self.db["structure_smiles"])

    def num_references(self):
        return len(self.db["reference_doi"])

    def num_couples(self):
        return len(self.db["c2t"])

    def num_couples_ref(self):
        return len(self.db["tc2r"])

    ### Taxonomy
    @functools.lru_cache(maxsize=None)
    def get_taxa(self) -> dict[int, str]:
        return self.db["taxonomy_names"]

    # TODO Not used
    def get_taxon_name_from_list_of_tids(self, tids: list[int]) -> list[str]:
        return [
            self.db["taxonomy_names"][tid]
            for tid in tids
            if tid in self.db["taxonomy_names"]
        ]

    def get_dict_of_tid_to_taxon_name(self, tid: Iterable[int]) -> dict[int, str]:
        return {
            t: self.db["taxonomy_names"][t]
            for t in tid
            if t in self.db["taxonomy_names"]
        }

    def get_taxon_name_from_tid(self, tid: int) -> str | None:
        try:
            tid = int(tid)
        except ValueError:
            return None
        if tid not in self.db["taxonomy_names"]:
            return None
        return self.db["taxonomy_names"][tid]

    def get_taxa_with_name_containing(self, query: str) -> Iterable[int]:
        query = query.lower()

        for tid, name in self.db["taxonomy_names"].items():
            if query in name.lower():
                yield tid

    # TODO Not used
    def get_rank_name_from_wid(self, wid: int) -> str | None:
        if wid not in self.db["taxonomy_ranks_names"]:
            return None
        return self.db["taxonomy_ranks_names"][wid]

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
        try:
            wid = int(wid)
        except ValueError:
            return None
        if wid in self.db["taxonomy_ranks"]:
            i_ranks = self.db["taxonomy_ranks"][wid]
            n_ranks = [self.get_rank_name_from_wid(int(it)) for it in i_ranks]
            ranks = " (" + ", ".join(set(n_ranks)) + ")"
        else:
            ranks = ""
        return ranks

    def get_taxonomic_tree(self, tid: int) -> list[tuple[int, int]]:
        if tid not in self.db["taxonomy_direct_parents"]:
            return []
        parent_taxa = self.db["taxonomy_direct_parents"][tid]
        tree = []
        for parent in parent_taxa:
            tree.append((parent, 1))
            if parent in self.db["taxonomy_parents_with_distance"]:
                for relative in self.db["taxonomy_parents_with_distance"][parent]:
                    distance = self.db["taxonomy_parents_with_distance"][parent][
                        relative
                    ]
                    tree.append((relative, distance))
        tree = sorted(tree, key=lambda x: x[1])
        return tree

    ### Structureonomy
    @functools.lru_cache(maxsize=None)
    def structures_set(self) -> set[int]:
        return set(self.db["structure_wid"])

    def get_structure_smiles_from_sid(self, sid: int) -> str | None:
        try:
            sid = self.db["structure_id"][sid]
            return self.db["structure_smiles"][sid]
        except (IndexError, ValueError):
            log.warning(f"Impossible to find a structure with sid={sid}")
            return None

    def get_structure_smiles_from_list_of_sids(self, sids: list[int]) -> list[str]:
        ids = [
            self.db["structure_id"][sid]
            for sid in sids
            if sid in self.db["structure_id"]
        ]
        llen = self.db["structure_smiles"]
        return [self.db["structure_smiles"][i] for i in ids if 0 <= i < len(llen)]

    def get_dict_of_sid_to_smiles(self, sid: Iterable[int]) -> dict[int, str]:
        ids = {
            s: self.db["structure_id"][s] for s in sid if s in self.db["structure_id"]
        }
        llen = self.db["structure_smiles"]
        return {
            sid: self.db["structure_smiles"][i]
            for sid, i in ids.items()
            if 0 <= i < len(llen)
        }

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

    def structure_get_tsv_from_scores(self, sids: list[int], scores) -> str:
        out = "Wikidata link\tSimilarity\tSmiles\n"
        for idx, score in enumerate(scores):
            sid = sids[idx]
            smiles = self.db["structure_smiles"][self.db["structure_id"][sid]]
            out += f"http://www.wikidata.org/entity/Q{sid}\t{score:.3f}\t{smiles}\n"
        return out

    ### Biblionomy
    @functools.lru_cache(maxsize=None)
    def get_refs(self) -> dict[int, str]:
        return self.db["reference_doi"]

    def get_ref_doi_from_list_of_rids(self, rids: list[int]) -> list[str]:
        return [
            self.db["reference_doi"][rid]
            for rid in rids
            if rid in self.db["reference_doi"]
        ]

    def get_dict_of_rid_to_ref_doi(self, rid: Iterable[int]) -> dict[int, str]:
        return {
            r: self.db["reference_doi"][r] for r in rid if r in self.db["reference_doi"]
        }

    def get_ref_doi_from_rid(self, rid: int) -> str | None:
        try:
            rid = int(rid)
        except ValueError:
            return None
        if rid not in self.db["reference_doi"]:
            return None
        return self.db["reference_doi"][rid]

    # TODO not sure it is the best way to proceed
    def get_references_with_doi(self, doi: str) -> Iterable[int]:
        for rid, d in self.db["reference_doi"].items():
            if doi == d:
                yield d

    ### Mixonomy
    def get_structures_of_taxon(self, tid: int, recursive: bool = True) -> list[int]:
        if tid in self.db["t2c"]:
            matching_structures = set(self.db["t2c"][tid])
        else:
            matching_structures = set()

        if recursive:
            if tid in self.db["taxonomy_children"]:
                for parent in self.db["taxonomy_children"][tid]:
                    if parent in self.db["t2c"]:
                        for structure in self.db["t2c"][parent]:
                            matching_structures.add(structure)

        return list(matching_structures)

    # TODO add get_numbers_of_structures_of_taxon

    def get_taxa_containing_structure(self, sid: int) -> set[int]:
        if sid in self.db["c2t"]:
            return self.db["c2t"][sid]
        return set()

    def get_number_of_taxa_containing_structure(self, sid: int) -> int:
        if sid not in self.db["c2t"]:
            return 0
        return len(self.db["c2t"][sid])

    ### WIP
    # def get_structures_of_reference(self, rid: int) -> list[int]:
    #     # TODO
    #     return list(matching_structures)

    # def get_taxa_of_reference(self, rid: int) -> list[int]:
    #     # TODO
    #     return list(matching_taxa)

    # def get_references_containing_couple(self, sid: int, tid: int) -> list[int]:
    #     # TODO
    #     return list(matching_couples)

    # def get_number_of_references_containing_couple(self, sid: int, tid: int) -> list[int]:
    #     # TODO
    #     return len(matching_couples)

    # def get_references_containing_structure(self, rid: int) -> list[int]:
    #     # TODO
    #     return list(matching_structures)

    # def get_number_of_references_containing_structure(self, rid: int) -> list[int]:
    #     # TODO
    #     return len(matching_structures)

    # def get_references_containing_taxa(self, rid: int) -> list[int]:
    #     # TODO
    #     return list(matching_taxa)

    # def get_number_of_references_containing_taxa(self, rid: int) -> list[int]:
    #     # TODO
    #     return len(matching_taxa)
