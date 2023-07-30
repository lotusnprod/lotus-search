# database.pkl

- **iid** (int)ernal compound id (can change at any update)
- **tid** (int): taxon id (wikidata id)
- **cid** (int): compound id (wikidata id)
- **rid** (int): reference id (wikidata id)
- **rank_id** (int): taxonomic rank id (wikidata id)


| key                            | type                            | content                                                           |
|--------------------------------|---------------------------------|-------------------------------------------------------------------|
| smileses                       | list[str]                       | Array of smiles strings                                           |
| sim_fps                        | list[bytes]                     | Similarity fingerprints                                           |
| sub_fps                        | list[bytes]                     | Substructure fingerprints                                         |
| t2c                            | dict[tid, set[cid]]             | All compounds from taxon tid                                      |
| c2t                            | dict[cid, set[tid]]             | All taxa containg compound cid                                    |
| tc2r                           | dict[tuple[tid, cid], set[rid]] | References of couple tid/cid                                      |
| compound_wid_to_id             | dict[cid, iid]                  | Translation of wids into iids                                     |
| compound_id_to_wid             | dict[iid, cid]                  | Translation of iids into cids                                     |
| taxonomy_direct_parents        | dict[tid, set[tid]]             | Direct taxonomic parents                                          |
| taxonomy_names                 | dict[tid, str]                  | Name of taxon tid                                                 |
| taxonomy_ranks                 | dict[tid, rank_id]              | Rank of taxon tid                                                 |
| taxonomy_children              | dict[tid, set[tid]]             | All the children of taxon tid                                     |
| taxonomy_parents_with_distance | dict[tid, dict[tid, int]]       | All the parents of taxon tid with an arbitrary taxonomic distance |
| taxonomy_ranks_names           | dict[tid, str]                  | Name of rank                                                      |
