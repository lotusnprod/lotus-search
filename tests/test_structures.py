from storage.models.structures import Structures


def test_structures():
    structure = Structures(id=1, smiles="C")
    assert (
        structure.__repr__()
        == "Structures(id=1, smiles=C, smiles_no_stereo=None, inchi=None, inchi_no_stereo=None, inchikey=None, inchikey_no_stereo=None, formula=None)"
    )
