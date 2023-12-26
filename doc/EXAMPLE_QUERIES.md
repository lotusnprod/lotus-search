## Structures

### Get all the structures that match the query

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/structures -d '{"structure_wid": 27151406, "taxon_wid": 158572 ,"reference_wid": 44488598}'
```

### Get all chlorinated structures

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/structures -d '{"structure":"Cl", "substructure_search": true, "taxon_name": "Gentiana septemfida", "reference_doi": "10.1021/NP50081A018"}'
```

## Taxa

### Get all the taxa that match the query

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/taxa -d '{"structure_wid": 27151406, "taxon_wid": 158572 ,"reference_wid": 44488598}'
```

### Get all the taxa where chlorinated structures are found in

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/taxa -d '{"structure":"Cl", "substructure_search": true, "taxon_name": "Gentiana septemfida", "reference_doi": "10.1021/NP50081A018"}'
```

## References

### Get all the references that match the query

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/references -d '{"structure_wid": 27151406, "taxon_wid": 158572 ,"reference_wid": 44488598}'
```

### Get all the references where chlorinated structures are found in

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/references -d '{"structure":"Cl", "substructure_search": true, "taxon_name": "Gentiana septemfida", "reference_doi": "10.1021/NP50081A018"}'
```

## Couples

### Get all the couples that match the query

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/couples -d '{"structure_wid": 27151406, "taxon_wid": 158572 ,"reference_wid": 44488598}'
```

### Get all the couples with chlorine

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/couples -d '{"structure":"Cl", "substructure_search": true, "taxon_name": "Gentiana septemfida", "reference_doi": "10.1021/NP50081A018"}'
```
