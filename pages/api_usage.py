import dash
import dash_bootstrap_components as dbc
from dash import dcc

# Register page (exposed in top navigation)
dash.register_page(
    __name__,
    path="/api",
    name="API usage",
    order=60,
    top_nav=True,
    title="LOTUS API usage",
)

_API_MD = r"""
# LOTUS Search API Usage

The REST API lets you programmatically query references, structures, taxa, triplets,
and utility endpoints (autocomplete, depiction, descriptors).

Base URL (default local dev): `http://localhost:8000`
(If using docker‑compose, API may be on `http://localhost:5000`.)

Versioning: endpoints are available at:
- `/v1/...` (stable explicit version)
- `/latest/...` (alias to most recent)
- Unversioned (e.g. `/structures`) also exists but use versioned paths for reliability.

## Endpoints Summary
| Method | Path | Description |
| ------ | ---- | ----------- |
| POST | `/v1/structures` | Search structures (similarity / substructure / formula / descriptors) |
| POST | `/v1/taxa` | Search taxa (by name or ID, optionally include children) |
| POST | `/v1/references` | Search references (by wid / doi / title / date range / journal) |
| POST | `/v1/triplets` | Intersect constraints and return (reference, structure, taxon) triplets |
| POST | `/v1/autocomplete/taxa` | Prefix autocomplete for taxon names |
| POST | `/v1/depiction/structure` | Return SVG depiction of a structure (optional highlight) |
| GET  | `/v1/descriptors/` | List available RDKit descriptor names |

All POST endpoints accept a JSON body shaped like the `Item` model. Unused
sections can be omitted.

## Core Request Model (abridged)
```json
{
  "reference": {
    "wid": null,
    "doi": null,
    "title": null,
    "option": {"date_min": null, "date_max": null, "journal": null}
  },
  "structure": {
    "wid": null,
    "molecule": null,
    "formula": null,
    "option": {
      "descriptors": null,
      "return_descriptors": false,
      "substructure_search": false,
      "similarity_level": 1.0,
      "sdf": false
    }
  },
  "taxon": {
    "wid": null,
    "name": null,
    "option": {"taxon_children": false}
  },
  "limit": null,
  "modeEnum": "objects"  // or "ids"
}
```

Notes:
- `limit = 0` means "no truncation" (return all matches).
- `modeEnum = "objects"` returns expanded objects; `"ids"` returns only identifiers.
- Provide **only one** of (wid / molecule / formula) for `structure`.
- Provide **only one** of (wid / doi / title) for `reference`.
- Provide **only one** of (wid / name) for `taxon`.

## Examples
### 1. Structure Similarity Search
```bash
curl -X POST http://localhost:8000/v1/structures \
  -H 'Content-Type: application/json' \
  -d '{"structure": {"molecule": "CCO", "option": {"substructure_search": false, "similarity_level": 0.7}}, "limit": 25, "modeEnum": "objects"}'
```

### 2. Substructure Search
```bash
curl -X POST http://localhost:8000/v1/structures \
  -H 'Content-Type: application/json' \
  -d '{"structure": {"molecule": "c1ccccc1", "option": {"substructure_search": true}}, "modeEnum": "ids"}'
```

### 3. Taxon Name Search (with children)
```bash
curl -X POST http://localhost:8000/v1/taxa \
  -H 'Content-Type: application/json' \
  -d '{"taxon": {"name": "Gentiana", "option": {"taxon_children": true}}, "modeEnum": "objects"}'
```

### 4. Reference Query by DOI (return objects)
```bash
curl -X POST http://localhost:8000/v1/references \
  -H 'Content-Type: application/json' \
  -d '{"reference": {"doi": "10.1080/1057563021000040466"}, "modeEnum": "objects"}'
```

### 5. Triplets Intersection (structure + taxon)
```bash
curl -X POST http://localhost:8000/v1/triplets \
  -H 'Content-Type: application/json' \
  -d '{"structure": {"molecule": "C", "option": {"substructure_search": true}}, "taxon": {"name": "lutea"}, "limit": 50, "modeEnum": "ids"}'
```

### 6. Autocomplete
```bash
curl -X POST http://localhost:8000/v1/autocomplete/taxa \
  -H 'Content-Type: application/json' \
  -d '{"taxon_name": "Gent"}'
```

### 7. Structure Depiction with Highlight
```bash
curl -X POST http://localhost:8000/v1/depiction/structure \
  -H 'Content-Type: application/json' \
  -d '{"structure": "CCOC(=O)C", "highlight": "COC"}'
```
Returns JSON with an inline SVG string.

### 8. Available Descriptors
```bash
curl http://localhost:8000/v1/descriptors/
```

## Error Responses
- 400: Invalid combination (e.g. both `wid` and `molecule`).
- 422: Invalid chemical string / formula / descriptor specification.
- 200 with empty lists: Valid query, no matches.

Example error:
```json
{"detail": "Only one of ['wid', 'molecule', 'formula'] should be provided"}
```

## Performance Tips
- Use `modeEnum = "ids"` for large exploratory queries, then fetch objects selectively.
- Narrow with date / journal / similarity filters early.
- Substructure queries can be slower; start broad with similarity first if possible.

## Integration (Python requests)
```python
import requests
payload = {"structure": {"molecule": "CCO"}, "modeEnum": "ids", "limit": 10}
r = requests.post("http://localhost:8000/v1/structures", json=payload, timeout=30)
r.raise_for_status()
print(r.json()["ids"])
```

## Stability & Backwards Compatibility
Versioned endpoints (`/v1/...`) will keep response field names and semantics stable.
Use `/latest/...` only if you are ready to pick up future versions automatically.

## License / Attribution
Data originates from LOTUS and Wikidata. Cite appropriately when publishing results.

---
Questions or suggestions? Open an issue on the repository.
"""

def layout():  # pragma: no cover - UI layer
    return dbc.Container([
        dcc.Markdown(_API_MD, link_target="_blank"),
    ])

