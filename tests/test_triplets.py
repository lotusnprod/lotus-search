from storage.models.triplets import Triplets


def test_triplets():
    triplets = Triplets(id=1, reference_id=1, structure_id=1, taxon_id=1)
    assert triplets.__repr__() == "Triplets(id=1, reference=1, structure=1, taxon=1)"
