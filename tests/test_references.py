from storage.references import References


def test_references():
    ref = References(id=1, doi="42.1/1")
    assert ref.__repr__() == "References(id=1, doi=42.1/1)"
