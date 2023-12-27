import sqlite3
from pathlib import Path
from typing import Any

from sqlalchemy import create_engine, text
from sqlalchemy.orm import sessionmaker

# Keep that this way so metadata gets all the tables
from storage.schemas import Base, SchemaVersion, Triplets


class Storage:
    SCHEMA_VERSION = 2

    def __init__(self, path: Path):
        self.con = sqlite3.connect(path / "index.db")
        self.list_limit = self.con.getlimit(sqlite3.SQLITE_LIMIT_VARIABLE_NUMBER)
        self.con.close()
        self.db_path = str((path / 'index.db').absolute())
        self.engine = create_engine(f"sqlite:///{self.db_path}")
        # Check if the schema table exists if not, call create
        new_db = self.query("SELECT name FROM sqlite_master WHERE type='table' AND name='schema_version'")

        if len(new_db) == 0:
            Base.metadata.create_all(self.engine)
            with self.session() as session:
                session.add(SchemaVersion(version=self.SCHEMA_VERSION))
                session.commit()

        with self.session() as session:
            db_version = session.query(SchemaVersion).first().version

        if db_version != self.SCHEMA_VERSION:
            raise Exception(
                f"Database schema version {db_version} does not match expected version {self.SCHEMA_VERSION}"
            )

    def session(self):
        Session = sessionmaker(bind=self.engine)
        return Session()

    def query(self, query: str):
        with self.engine.connect() as connection:
            result = connection.execute(text(query))
        return result.fetchall()

    def add_triplets(self, triplets: list[list]) -> None:
        with self.session() as session:
            mappings = [{"reference_id": t[0], "structure_id": t[1], "taxon_id": t[2]} for t in triplets]
            session.bulk_insert_mappings(Triplets, mappings)
            session.commit()

    def get_generic_of_generic(self, out: Any, inp: Any, item: int) -> set[int]:
        with self.session() as session:
            result = session.query(out).filter(inp == item).distinct()
            return {row[0] for row in result}

    def get_references_containing_taxon(self, tid: int) -> set[int]:
        return self.get_generic_of_generic(Triplets.reference_id, Triplets.taxon_id, tid)

    def get_references_containing_structure(self, sid: int) -> set[int]:
        return self.get_generic_of_generic(Triplets.reference_id, Triplets.structure_id, sid)

    def get_taxa_of_reference(self, rid: int) -> set[int]:
        return self.get_generic_of_generic(Triplets.taxon_id, Triplets.reference_id, rid)

    def get_structures_of_reference(self, rid: int) -> set[int]:
        return self.get_generic_of_generic(Triplets.structure_id, Triplets.reference_id, rid)

    def get_taxa_of_structure(self, sid: int) -> set[int]:
        return self.get_generic_of_generic(Triplets.taxon_id, Triplets.structure_id, sid)

    def get_structures_of_taxon(self, tid: int) -> set[int]:
        return self.get_generic_of_generic(Triplets.structure_id, Triplets.taxon_id, tid)

    def validate_column_name(self, s):
        ps = s.replace('_', '')
        return ps.islower() and ps.isalpha()

    def get_generics_of_generics(self, out: Any, inp: Any, items: set[int]) -> set[int]:
        if len(items) > self.list_limit:
            raise Exception(
                f"Too many {inp}, the limit is {self.list_limit} and you have {len(items)}"
            )

        with self.session() as session:
            result = session.query(out).filter(inp.in_(items)).distinct().all()
            return {row[0] for row in result}

    def get_references_of_structures(self, structures: set[int]) -> set[int]:
        return self.get_generics_of_generics(Triplets.reference_id, Triplets.structure_id, structures)

    def get_references_of_taxa(self, taxa: set[int]) -> set[int]:
        return self.get_generics_of_generics(Triplets.reference_id, Triplets.taxon_id, taxa)

    def get_structures_of_references(self, references: set[int]) -> set[int]:
        return self.get_generics_of_generics(Triplets.structure_id, Triplets.reference_id, references)

    def get_taxa_of_structures(self, structures: set[int]) -> set[int]:
        return self.get_generics_of_generics(Triplets.taxon_id, Triplets.structure_id, structures)

    def get_taxa_of_references(self, references: set[int]) -> set[int]:
        return self.get_generics_of_generics(Triplets.taxon_id, Triplets.reference_id, references)

    def get_triplets_for(self,
                         reference_ids: set[int] | None,
                         structure_ids: set[int] | None,
                         taxon_ids: set[int] | None) -> set[tuple[int, int, int]]:
        with self.session() as session:
            filters = []
            if reference_ids is not None:
                filters += [Triplets.reference_id.in_(reference_ids)]
            if structure_ids is not None:
                filters += [Triplets.structure_id.in_(structure_ids)]
            if taxon_ids is not None:
                filters += [Triplets.taxon_id.in_(taxon_ids)]
            result = session.query(Triplets.reference_id,
                                   Triplets.structure_id,
                                   Triplets.taxon_id)

            if len(filters) > 0:
                result = result.filter(*filters)

            return {(row[0], row[1], row[2]) for row in result.distinct()}
