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

    def get_taxon_name_from_list_of_wid(self, wid: list[int]) -> list[str]:
        return [
            self.db["taxonomy_names"][w] for w in wid if w in self.db["taxonomy_names"]
        ]

    def get_dict_of_wid_to_taxon_name(self, wid: Iterable[int]) -> dict[int, str]:
        return {
            w: self.db["taxonomy_names"][w]
            for w in wid
            if w in self.db["taxonomy_names"]
        }

    def get_taxon_name_from_wid(self, wid: int) -> str | None:
        try:
            wid = int(wid)
        except ValueError:
            return None
        if wid not in self.db["taxonomy_names"]:
            return None
        return self.db["taxonomy_names"][wid]

    def get_taxa_with_name_containing(self, query: str) -> Iterable[int]:
        query = query.lower()

        for wid, name in self.db["taxonomy_names"].items():
            if query in name.lower():
                yield wid

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

    def get_taxonomic_tree(self, wid: int) -> list[tuple[int, int]]:
        if wid not in self.db["taxonomy_direct_parents"]:
            return []
        parent_taxa = self.db["taxonomy_direct_parents"][wid]
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

    def get_structure_smiles_from_wid(self, wid: int) -> str | None:
        try:
            sid = self.db["structure_id"][wid]
            return self.db["structure_smiles"][sid]
        except (IndexError, ValueError):
            log.warning(f"Impossible to find a structure with wid={wid}")
            return None

    def get_structure_smiles_from_list_of_wid(self, wid: list[int]) -> list[str]:
        ids = [self.db["structure_id"][w] for w in wid if w in self.db["structure_id"]]
        llen = self.db["structure_smiles"]
        return [self.db["structure_smiles"][i] for i in ids if 0 <= i < len(llen)]

    def get_dict_of_wid_to_smiles(self, wid: Iterable[int]) -> dict[int, str]:
        ids = {
            w: self.db["structure_id"][w] for w in wid if w in self.db["structure_id"]
        }
        llen = self.db["structure_smiles"]
        return {
            wid: self.db["structure_smiles"][i]
            for wid, i in ids.items()
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

        new_keys = [self.db["structure_wid"][iid] for iid in iids]
        out = []
        for iid, wid in zip(iids, new_keys):
            out.append((wid, DataStructs.TanimotoSimilarity(fp, fp_db[iid])))
        return out

    def structure_get_tsv_from_scores(self, wids, scores) -> str:
        out = "Wikidata link\tSimilarity\tSmiles\n"
        for idx, score in enumerate(scores):
            wid = wids[idx]
            smiles = self.db["structure_smiles"][self.db["structure_id"][wid]]
            out += f"http://www.wikidata.org/entity/Q{wid}\t{score:.3f}\t{smiles}\n"
        return out

    ### Biblionomy
    @functools.lru_cache(maxsize=None)
    def get_refs(self) -> dict[int, str]:
        return self.db["reference_doi"]

    def get_ref_doi_from_list_of_wid(self, wid: list[int]) -> list[str]:
        return [
            self.db["reference_doi"][w] for w in wid if w in self.db["reference_doi"]
        ]

    def get_dict_of_wid_to_ref_doi(self, wid: Iterable[int]) -> dict[int, str]:
        return {
            w: self.db["reference_doi"][w] for w in wid if w in self.db["reference_doi"]
        }

    def get_ref_doi_from_wid(self, wid: int) -> str | None:
        try:
            wid = int(wid)
        except ValueError:
            return None
        if wid not in self.db["reference_doi"]:
            return None
        return self.db["reference_doi"][wid]

    # TODO not sure it is the best way to proceed
    def get_references_with_doi(self, doi: str) -> Iterable[int]:
        for wid, d in self.db["reference_doi"].items():
            if doi == d:
                yield d

    ### Mixonomy
    def get_structures_of_taxon(self, wid: int, recursive: bool = True) -> list[int]:
        if wid in self.db["t2c"]:
            matching_structures = set(self.db["t2c"][wid])
        else:
            matching_structures = set()

        if recursive:
            if wid in self.db["taxonomy_children"]:
                for parent in self.db["taxonomy_children"][wid]:
                    if parent in self.db["t2c"]:
                        for structure in self.db["t2c"][parent]:
                            matching_structures.add(structure)

        return list(matching_structures)

    # TODO add get_numbers_of_structures_of_taxon

    def get_taxa_containing_structure(self, wid: int) -> set[int]:
        if wid in self.db["c2t"]:
            return self.db["c2t"][wid]
        return set()

    def get_number_of_taxa_containing_structure(self, wid: int) -> int:
        if wid not in self.db["c2t"]:
            return 0
        return len(self.db["c2t"][wid])

    ### WIP
    # def get_structures_of_reference(self, wid: int) -> list[int]:
    #     # TODO
    #     return list(matching_structures)

    # def get_taxa_of_reference(self, wid: int) -> list[int]:
    #     # TODO
    #     return list(matching_taxa)

    # def get_references_containing_couple(self, sid: int, tid: int) -> list[int]:
    #     # TODO
    #     return list(matching_couples)

    # def get_number_of_references_containing_couple(self, sid: int, tid: int) -> list[int]:
    #     # TODO
    #     return len(matching_couples)

    # def get_references_containing_structure(self, wid: int) -> list[int]:
    #     # TODO
    #     return list(matching_structures)

    # def get_number_of_references_containing_structure(self, wid: int) -> list[int]:
    #     # TODO
    #     return len(matching_structures)

    # def get_references_containing_taxa(self, wid: int) -> list[int]:
    #     # TODO
    #     return list(matching_taxa)

    # def get_number_of_references_containing_taxa(self, wid: int) -> list[int]:
    #     # TODO
    #     return len(matching_taxa)
