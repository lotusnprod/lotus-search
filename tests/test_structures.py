from storage.structures import Structures


def test_structures():
    structure = Structures(id=1, smiles="C")
    assert structure.__repr__() == "Structures(id=1, smiles=C)"
