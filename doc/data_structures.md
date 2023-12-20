## database.pkl

- **iid** (int)ernal structure id (can change at any update)
- **tid** (int): taxon id (wikidata id)
- **sid** (int): structure id (wikidata id)
- **rid** (int): reference id (wikidata id)
- **rank_id** (int): taxonomic rank id (wikidata id)

| key                            | type                            | content                                                           |
|--------------------------------|---------------------------------|-------------------------------------------------------------------|
| structure_smiles               | list[str]                       | Array of smiles strings                                           |
| structure_wid                  | list[int]                       | Array of sids                                                     |
| structure_sim_fps              | list[bytes]                     | Similarity fingerprints                                           |
| structure_library              | bytes                           | The RDKit library                                                 |
| structure_sim_h_fps            | list[bytes]                     | Similarity fingerprints (with explicit Hs)                        |
| structure_library_h            | bytes                           | The RDKit library (with explicit Hs)                              |
| t2s                            | dict[tid, set[sid]]             | All structures from taxon tid                                     |
| s2t                            | dict[sid, set[tid]]             | All taxa containg structure sid                                   |
| ts2r                           | dict[tuple[tid, sid], set[rid]] | References of couple tid/sid                                      |
| taxonomy_direct_parents        | dict[tid, set[tid]]             | Direct taxonomic parents                                          |
| taxonomy_names                 | dict[tid, str]                  | Name of taxon tid                                                 |
| taxonomy_ranks                 | dict[tid, set[rank_id]]         | Rank of taxon tid                                                 |
| taxonomy_children              | dict[tid, set[tid]]             | All the children of taxon tid                                     |
| taxonomy_parents_with_distance | dict[tid, dict[tid, int]]       | All the parents of taxon tid with an arbitrary taxonomic distance |
| taxonomy_ranks_names           | dict[tid, str]                  | Name of rank                                                      |
| reference_doi                  | dict[rid, str]                  | Reference DOI                                                     |
