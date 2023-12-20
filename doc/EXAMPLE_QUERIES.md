## Structures

### Get all the structures that match the query

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/structures -d '{"structure":"C=C[C@@H]1[C@@H]2CCOC(=O)C2=CO[C@H]1O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)O", "taxon_name":"Gentiana lutea"}'
```

### Get all chlorinated structures

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/structures -d '{"structure":"Cl", "substructure_search": true}'
```

## Taxa

### Get all the taxa that match the query

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/taxa -d '{"structure":"C=C[C@@H]1[C@@H]2CCOC(=O)C2=CO[C@H]1O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)O", "taxon_name":"Gentiana lutea"}'
```

### Get all the taxa where chlorinated structures are found in

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/taxa -d '{"structure":"Cl", "substructure_search": true}'
```

## References

### Get all the references that match the query

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/references -d '{"structure":"C=C[C@@H]1[C@@H]2CCOC(=O)C2=CO[C@H]1O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)O", "taxon_name":"Gentiana lutea"}'
```

### Get all the references where chlorinated structures are found in

### TODO This one is taking super long, don't

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/references -d '{"structure":"Cl", "substructure_search": true}'
```

## Couples

### Get all the couples that match the query

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/couples -d '{"structure":"C=C[C@@H]1[C@@H]2CCOC(=O)C2=CO[C@H]1O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)O", "taxon_name":"Gentiana lutea"}'
```

### Get all the couples with chlorine

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/couples -d '{"structure":"Cl", "substructure_search": true}'
```