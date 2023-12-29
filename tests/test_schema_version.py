from storage.schema_version import SchemaVersion


def test_schema_version():
    schema_version = SchemaVersion(version=1)
    assert schema_version.__repr__() == "SchemaVersion(version=1)"
