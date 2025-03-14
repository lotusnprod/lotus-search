## Description

This application is available on [Nprod.net](https://search.nprod.net)

The dataset is from the [LOTUS](https://lotus.nprod.net) initiative and [Wikidata](https://www.wikidata.org).

## Requirements

- [Dash](https://dash.plotly.com) for its web ui.
- [EPAM's Ketcher](https://lifescience.opensource.epam.com/ketcher/index.html?ref=search.nprod.net) for the molecule
  editor.
- [Rdkit](https://www.rdkit.org) for the molecule massaging.

## Install & use

### The docker compose way (recommended)

```shell
docker-compose run -it backend python update.py  # You do not have to run it everytime, but that's recommended until we have a schema versioning
docker-compose up --build
```

The web server is then available on <http://localhost:3000> and the API on <http://localhost:5000>.

### Run tests
```
make tests
```

### The manual way

To run it yourself, the source is available at: <https://github.com/lotusnprod/lotus-search>:

- Install dependencies using uv
    - If you do not have uv installed:
        - `curl -LsSf https://astral.sh/uv/install.sh | sh`
    - Then:
        - `uv sync`
- Run `python update.py` (takes a few minutes)
- Run `uvicorn main:app --reload` (almost instant)

## Authors

[Adriano Rutz](https://adafede.github.io), is the one that pushed me into that, we are both part of the team behind
LOTUS.

### Jonathan Bisson

- Personal website: <https://www.bjonnh.net>
- You can find me on mastodon too(t): <https://mastodon.social/@bjonnh>

### Adriano Rutz

- Personal website: <https://adafede.github.io>
- And on mastodon: <https://mastodon.online/@adafede>

## Data safety

Your molecules are never stored unless they make our service crash.

But if your molecules are super secret, it is like your extremities, don't insert
them in machines you don't know or understand.

Also, this is an experimental tool meant to test things,
you're not supposed to rely on it for anything important, and
it comes with no warranty or support whatsoever.

## License and legalese

<https://raw.githubusercontent.com/lotusnprod/lotus-search/main/LICENSE>
