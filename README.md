This application is available on [Nprod.net](https://search.nprod.net/)

The code of this is at: https://github.com/lotusnprod/text_search

The dataset is from the [LOTUS](https://lotus.nprod.net/) database and [Wikidata](https://www.wikidata.org).

[Adriano Rutz](https://adafede.github.io/), is the one that pushed me into that, we are both part of the team behind LOTUS.

It is using [Streamlit](https://streamlit.io)  for its web ui.
And [EPAM's Ketcher](https://lifescience.opensource.epam.com/ketcher/index.html?ref=blog.streamlit.io) for the molecule editor.

[Rdkit](https://www.rdkit.org) is used for the molecule massaging.

To run it yourself, the source is on [GitHub](https://github.com/lotusnprod/mol_search):
- Install and run poetry update
- Run `get_smiles_wikidata.py`  (takes 20s)
- Run `generate_database.py`    (takes 1min on a Ryzen 7 3700X)
- Run `streamlit run main.py`   (almost instant)

## **Me**

- Personal website: https://www.bjonnh.net
- You can find me on mastodon too(t): https://mastodon.social/@bjonnh

## **Data safety**

Your molecules are never stored, the only people that could eventually see them are the streamlit people. 
And I doubt they care about your molecules. But if your molecules are super secret, it is like your extremities, don't insert 
them in machines you don't know or understand.

Also, this is an experimental tool meant to test things, you're not supposed to rely on it for anything important, and
it comes with no warranty or support whatsoever.

## **License and legalese**

https://raw.githubusercontent.com/lotusnprod/text_search/main/LICENSE
