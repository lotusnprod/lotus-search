import sqlite3
from pathlib import Path


class Storage:
    SCHEMA_VERSION = 1

    def __init__(self, path: Path):
        self.con = sqlite3.connect(path / "index.db")
        self.list_limit = self.con.getlimit(sqlite3.SQLITE_LIMIT_VARIABLE_NUMBER)
        cursor = self.cursor()
        # Check if the schema table exists if not, call create
        cursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='schema'"
        )
        new_db = cursor.fetchone()
        if new_db is None:
            self.create_schema()

        cursor.execute("SELECT version FROM schema")
        db_version = cursor.fetchone()[0]
        if db_version != self.SCHEMA_VERSION:
            raise Exception(
                f"Database schema version {db_version} does not match expected version {self.SCHEMA_VERSION}"
            )

    def cursor(self):
        return self.con.cursor()

    def close(self):
        self.con.close()

    def create_schema(self):
        # Create the schema table, it contains only a value which is the current version

        self.create_version_table()

        self.create_triplets_table()

        self.con.commit()

    def create_version_table(self):
        cursor = self.cursor()
        cursor.execute("CREATE TABLE schema (version INTEGER)")
        cursor.execute("INSERT INTO schema VALUES (?)", (self.SCHEMA_VERSION,))
        self.con.commit()

    def create_triplets_table(self):
        cursor = self.cursor()
        cursor.execute(
            "CREATE TABLE triplets (id INTEGER PRIMARY KEY, reference_id INTEGER, structure_id INTEGER, taxon_id INTEGER)"
        )
        cursor.execute("CREATE INDEX triplets_reference_id ON triplets (reference_id)")
        cursor.execute("CREATE INDEX triplets_structure_id ON triplets (structure_id)")
        cursor.execute("CREATE INDEX triplets_taxon_id ON triplets (taxon_id)")
        cursor.execute(
            "CREATE INDEX triplets_reference_structure_id ON triplets (reference_id, structure_id)"
        )
        cursor.execute(
            "CREATE INDEX triplets_reference_taxon_id ON triplets (reference_id, taxon_id)"
        )
        cursor.execute(
            "CREATE INDEX triplets_structure_taxon_id ON triplets (structure_id, taxon_id)"
        )
        cursor.execute(
            "CREATE INDEX triplets_reference_structure_taxon_id ON triplets (reference_id, structure_id, taxon_id)"
        )
        self.con.commit()

    def add_triplets(self, triplets: list[list]) -> None:
        cursor = self.cursor()
        cursor.executemany(
            "INSERT INTO triplets (reference_id, structure_id, taxon_id) VALUES (?, ?, ?)",
            triplets,
        )
        self.con.commit()

    def get_references_containing_taxon(self, tid: int) -> set[int]:
        cursor = self.cursor()
        cursor.execute(
            "SELECT DISTINCT reference_id FROM triplets WHERE taxon_id = ?", (tid,)
        )
        return {row[0] for row in cursor.fetchall()}

    def get_references_containing_structure(self, sid: int) -> set[int]:
        cursor = self.cursor()
        cursor.execute(
            "SELECT DISTINCT reference_id FROM triplets WHERE structure_id = ?", (sid,)
        )
        return {row[0] for row in cursor.fetchall()}

    def get_references_of_couple(self, sid: int, tid: int) -> set[int]:
        cursor = self.cursor()
        cursor.execute(
            "SELECT DISTINCT reference_id FROM triplets WHERE structure_id = ? AND taxon_id = ?",
            (sid, tid),
        )
        return {row[0] for row in cursor.fetchall()}

    def get_taxa_of_reference(self, rid: int) -> set[int]:
        cursor = self.cursor()
        cursor.execute("SELECT DISTINCT taxon_id FROM triplets WHERE reference_id = ?", (rid,))
        return {row[0] for row in cursor.fetchall()}

    def get_structures_of_reference(self, rid: int) -> set[int]:
        cursor = self.cursor()
        cursor.execute(
            "SELECT DISTINCT structure_id FROM triplets WHERE reference_id = ?", (rid,)
        )
        return {row[0] for row in cursor.fetchall()}

    def get_taxa_of_structure(self, sid: int) -> set[int]:
        cursor = self.cursor()
        cursor.execute("SELECT DISTINCT taxon_id FROM triplets WHERE structure_id = ?", (sid,))
        return {row[0] for row in cursor.fetchall()}

    def get_structures_of_taxon(self, tid: int) -> set[int]:
        cursor = self.cursor()
        cursor.execute("SELECT DISTINCT structure_id FROM triplets WHERE taxon_id = ?", (tid,))
        return {row[0] for row in cursor.fetchall()}

    def validate_column_name(self, s):
        ps = s.replace('_', '')
        return ps.islower() and ps.isalpha()

    def get_generics_of_generics(self, inp: str, out: str, items: set[int]) -> set[int]:
        # You need to make absolutely sure that inp and out are not given by users
        # Check if inp and out contain only lowercase chars and underscore
        if not self.validate_column_name(inp) or not self.validate_column_name(out):
            raise Exception("Invalid input or output column name: in: {inp} out: {out}")

        # TODO less generic exception and catch it in the API
        # TODO better yet, iterate over chunks of list_limit and combine
        if len(items) > self.list_limit:
            raise Exception(
                f"Too many {inp}, the limit is {self.list_limit} and you have {len(items)}"
            )

        cursor = self.cursor()
        cursor.execute(
            f"SELECT DISTINCT {out} FROM triplets WHERE {inp} IN ({{}})".format(
                ",".join("?" * len(items))
            ),
            tuple(items),
        )

        return {row[0] for row in cursor.fetchall()}

    def get_references_of_structures(self, structures: set[int]) -> set[int]:
        return self.get_generics_of_generics("structure_id", "reference_id", structures)

    def get_references_of_taxa(self, taxa: set[int]) -> set[int]:
        return self.get_generics_of_generics("taxon_id", "reference_id", taxa)

    def get_structures_of_references(self, references: set[int]) -> set[int]:
        return self.get_generics_of_generics("reference_id", "structure_id", references)

    def get_taxa_of_structures(self, structures: set[int]) -> set[int]:
        return self.get_generics_of_generics("structure_id", "taxon_id", structures)

    def get_taxa_of_references(self, references: set[int]) -> set[int]:
        return self.get_generics_of_generics("reference_id", "taxon_id", references)
