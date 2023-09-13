## Structures
### Get all the molecules that contain a chlorine from every specie with taxus in the name 

Careful that's taxus anywhere in the name, so you'll get Cephalotaxus as well!

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:8000/v1_0/structures -d '{"molecule":"Cl", "substructure_search": true, "taxon_name":"Taxus"}'
```
## Taxa
### Get all the taxa with gent in their name producing a compound with chlorine in it

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:8000/v1_0/taxa -d '{"molecule":"Cl", "substructure_search": true, "taxon_name":"gent"}'
```

### Get all the taxa producing chlorinated compounds
```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:8000/v1_0/taxa -d '{"molecule":"Cl", "substructure_search": true}'
```
## Couples
### Get all the couples that contain that specific compound and that have taxus in the organism name
Careful that's taxus anywhere in the name, so you'll get Cephalotaxus as well!

```shell
curl -XPOST -H 'Accept: application/json' -H 'Content-Type: application/json' http://127.0.0.1:8000/v1_0/couples -d '{"molecule":"O=C(OC1C(=C)C2CC3(O)CC(OC(=O)C)C(=C(C(OC(=O)C)C(OC(=O)C)C2(C)C(OC(=O)C)C1)C3(C)C)C)C=CC=4C=CC=CC4", "taxon_name":"tax"}'
```
