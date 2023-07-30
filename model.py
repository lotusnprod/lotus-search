from collections.abc import Iterable
from typing import Any

from rdkit import Chem, DataStructs


class DataModel:
    def __init__(self, db):
        self.db = db

    def num_taxa(self):
        return len(self.db["taxonomy_names"])

    def num_compounds(self):
        return len(self.db["smileses"])

    def num_couples(self):
        return len(self.db["c2t"])

    ### Taxonomy
    def get_taxon_name_from_wid(self, wid: int) -> str | None:
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

    ### Compoundonomy
    def get_compound_wid_from_id(self, iid: int) -> int | None:
        try:
            return self.db["compound_id_to_wid"][iid]
        except ValueError:
            return None

    def get_compound_smiles_from_iid(self, iid: int) -> str | None:
        try:
            return self.db["smileses"][iid]
        except:
            print("Impossible to find a compound with wid={wid} cid={cid}")
            return None

    def get_compound_smiles_from_wid(self, wid: int) -> str | None:
        if wid not in self.db["compound_wid_to_id"]:
            return None
        cid = self.db["compound_wid_to_id"][wid]
        return self.get_compound_smiles_from_iid(cid)

    def compound_search(self, fp: bytes) -> list[tuple[int, float]]:
        results = DataStructs.BulkTanimotoSimilarity(fp, self.db["sim_fps"])
        return [(j, score) for j, score in enumerate(results)]

    def compound_search_substructure(self, fp: bytes, mol: Any) -> list[tuple[int, float]]:
        out = []
        for j in self.db["library"].GetMatches(mol, numThreads=-1, maxResults=-1):
            out.append((j, DataStructs.TanimotoSimilarity(fp, self.db["sim_fps"][j])))
        return out

    def compound_get_tsv_from_scores(self, scores) -> str:
        out = "Wikidata link\tSimilarity\tSmiles\n"
        for score in scores:
            wid = self.db["compound_id_to_wid"][score[0]]
            smiles = self.db["smileses"][score[0]]
            out += f"http://www.wikidata.org/entity/Q{wid}\t{score[1]}\t{smiles}\n"
        return out

    ### Taxonomy to compoundonomy
    def get_compounds_of_taxon(self, wid: int, recursive: bool = True) -> set[int]:
        if wid in self.db["t2c"]:
            matching_compounds = set(self.db["t2c"][wid])
        else:
            matching_compounds = set()

        if recursive:
            if wid in self.db["taxonomy_children"]:
                for parent in self.db["taxonomy_children"][wid]:
                    if parent in self.db["t2c"]:
                        for compound in self.db["t2c"][parent]:
                            matching_compounds.add(compound)

        return matching_compounds

    def get_taxa_containing_compound(self, wid: int) -> set[int]:
        if wid in self.db["c2t"]:
            return self.db["c2t"][wid]
        return set()

    def get_number_of_taxa_containing_compound(self, wid: int) -> int:
        if wid not in self.db["c2t"]:
            return 0
        return len(self.db["c2t"][wid])

    def get_ranks_string(self, wid: int) -> str:
        if wid in self.db["taxonomy_ranks"]:
            i_ranks = self.db["taxonomy_ranks"][wid]
            n_ranks = [self.get_rank_name_from_wid(int(it)) for it in i_ranks]
            ranks = " (" + ", ".join(set(n_ranks)) + ")"
        else:
            ranks = ""
        return ranks
