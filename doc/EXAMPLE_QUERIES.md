# Example queries

## Structures

### Get all the structures that match the query

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/structures -d '{"structure": {"wid": 27151406}, "taxon": {"wid": 158572}, "reference": {"wid":44488598}}'
```

### Get all chlorinated structures

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/structures -d '{"structure": {"molecule": "Cl", "option": {"substructure_search": true}}, "taxon": {"name": "Gentiana septemfida"}, "reference": {"doi": "10.1021/NP50081A018"}}'
```

## Taxa

### Get all the taxa that match the query

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/taxa -d '{"structure": {"wid": 27151406}, "taxon": {"wid": 158572}, "reference": {"wid":44488598}}'
```

### Get all the taxa where chlorinated structures are found in

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/taxa -d '{"structure": {"molecule": "Cl", "option": {"substructure_search": true}}, "taxon": {"name": "Gentiana septemfida"}, "reference": {"doi": "10.1021/NP50081A018"}}'
```

## References

### Get all the references that match the query

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/references -d '{"structure": {"wid": 27151406}, "taxon": {"wid": 158572}, "reference": {"wid":44488598}}'
```

### Get all the references where chlorinated structures are found in

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/references -d '{"structure": {"molecule": "Cl", "option": {"substructure_search": true}}, "taxon": {"name": "Gentiana septemfida"}, "reference": {"doi": "10.1021/NP50081A018"}}'
```

## Triplets

### Get all the triplets that match the query

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/triplets -d '{"structure": {"wid": 27151406}, "taxon": {"wid": 158572}, "reference": {"wid":44488598}}'
```

### Get all the triplets with chlorine

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:5000/v1_0/triplets -d '{"structure": {"molecule": "Cl", "option": {"substructure_search": true}}, "taxon": {"name": "Gentiana septemfida"}, "reference": {"doi": "10.1021/NP50081A018"}}'
```

### TODO add new ones to showcase formula/descriptors/journal search
