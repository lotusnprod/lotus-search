from collections.abc import Iterable
from typing import Any

from rdkit import Chem, DataStructs


class DataModel:
    def __init__(self, db):
        self.db = db

    def num_taxa(self):
        return len(self.db["taxonomy_names"])

    def num_compounds(self):
        return len(self.db["compound_smiles"])

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
    def get_compound_smiles_from_wid(self, wid: int) -> str | None:

        try:
            cid = self.db["compound_wid"].index(wid)
            return self.db["compound_smiles"][cid]
        except (IndexError, ValueError):
            print(f"Impossible to find a compound with wid={wid}")
            return None

    def get_compound_smiles_from_list_of_wid(self, wid: list[int]) -> list[str]:
        ids = [self.db["compound_id"][w] for w in wid if w in self.db["compound_id"]]
        llen = self.db["compound_smiles"]
        return [self.db["compound_smiles"][i] for i in ids if i >= 0 and i < len(llen)]

    def compound_search(self, fp: bytes) -> list[tuple[int, float]]:
        scores = DataStructs.BulkTanimotoSimilarity(fp, self.db["compound_sim_fps"])
        return [(wid, score) for wid, score in zip(self.db["compound_wid"], scores)]

    def compound_search_substructure(self, fp: bytes, mol: Any, chirality: bool) -> list[tuple[int, float]]:
        out = []
        iids = self.db["compound_library"].GetMatches(mol, numThreads=-1, maxResults=-1, useQueryQueryMatches=True,
                                                      useChirality=chirality)

        new_keys = [self.db["compound_wid"][iid] for iid in iids]
        for iid, wid in zip(iids, new_keys):
            out.append((wid, DataStructs.TanimotoSimilarity(fp, self.db["compound_sim_fps"][iid])))
        return out

    def compound_get_tsv_from_scores(self, scores) -> str:
        out = "Wikidata link\tSimilarity\tSmiles\n"
        iids = {score[0]: self.db["compound_wid"].index(score[0]) for score in scores}
        for score in scores:
            wid = score[0]
            smiles = self.db["compound_smiles"][iids[score[0]]]
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
