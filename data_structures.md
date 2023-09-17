# database.pkl

- **iid** (int)ernal compound id (can change at any update)
- **tid** (int): taxon id (wikidata id)
- **cid** (int): compound id (wikidata id)
- **rid** (int): reference id (wikidata id)
- **rank_id** (int): taxonomic rank id (wikidata id)


| key                            | type                            | content                                                           |
|--------------------------------|---------------------------------|-------------------------------------------------------------------|
| compound_smiles                | list[str]                       | Array of smiles strings                                           |
| compound_wid                   | list[int]                       | Array of cids                                                     |
| compound_sim_fps               | list[bytes]                     | Similarity fingerprints                                           |
| compound_library               | bytes                           | The RDKit library                                                 |
| compound_sim_h_fps             | list[bytes]                     | Similarity fingerprints (with explicit Hs)                        |
| compound_library_h             | bytes                           | The RDKit library (with explicit Hs)                              |
| t2c                            | dict[tid, set[cid]]             | All compounds from taxon tid                                      |
| c2t                            | dict[cid, set[tid]]             | All taxa containg compound cid                                    |
| tc2r                           | dict[tuple[tid, cid], set[rid]] | References of couple tid/cid                                      |
| taxonomy_direct_parents        | dict[tid, set[tid]]             | Direct taxonomic parents                                          |
| taxonomy_names                 | dict[tid, str]                  | Name of taxon tid                                                 |
| taxonomy_ranks                 | dict[tid, rank_id]              | Rank of taxon tid                                                 |
| taxonomy_children              | dict[tid, set[tid]]             | All the children of taxon tid                                     |
| taxonomy_parents_with_distance | dict[tid, dict[tid, int]]       | All the parents of taxon tid with an arbitrary taxonomic distance |
| taxonomy_ranks_names           | dict[tid, str]                  | Name of rank                                                      |
