from storage.models.references import References


def test_references():
    ref = References(
        id=1,
        doi="42.1/1",
        title="Iridoids from Seeds of Gentiana Lutea",
        date="2003-08",
        journal="Natural Product Research",
    )
    assert (
        ref.__repr__()
        == "References(id=1, doi=42.1/1, title=Iridoids from Seeds of Gentiana Lutea, date=2003-08, journal=Natural Product Research)"
    )
