from storage.models.journals import Journals


def test_journals():
    ref = Journals(
        id=1,
        title="Journal A",
    )
    assert ref.__repr__() == "Journals(id=1, title=Journal A)"
