import sqlite3
from pathlib import Path
from typing import Any

from sqlalchemy import create_engine, text
from sqlalchemy.dialects.sqlite import insert as sqlite_upsert
from sqlalchemy.orm import sessionmaker

# Keep that this way so metadata gets all the tables
from storage.models import (
    Base,
    References,
    SchemaVersion,
    Structures,
    TaxoNames,
    TaxoParents,
    TaxoRankNames,
    TaxoRanks,
    Triplets,
)


class Storage:
    SCHEMA_VERSION = 4

    def __init__(self, path: Path):
        self.con = sqlite3.connect(path / "index.db")
        self.list_limit = self.con.getlimit(sqlite3.SQLITE_LIMIT_VARIABLE_NUMBER)
        self.con.close()
        self.db_path = str((path / "index.db").absolute())
        self.engine = create_engine(f"sqlite:///{self.db_path}")
        # Check if the schema table exists if not, call create
        new_db = self.query(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='schema_version'"
        )

        if len(new_db) == 0:
            Base.metadata.create_all(self.engine)
            with self.session() as session:
                session.add(SchemaVersion(version=self.SCHEMA_VERSION))
                session.commit()

        with self.session() as session:
            db_version = session.query(SchemaVersion).first().version

        if db_version != self.SCHEMA_VERSION:
            raise Exception(
                f"Database schema version {db_version} does not match expected version {self.SCHEMA_VERSION}. You want to delete data/index.db and rebuild it."
            )

    def session(self):
        Session = sessionmaker(bind=self.engine)
        return Session()

    def query(self, query: str):
        with self.engine.connect() as connection:
            result = connection.execute(text(query))
        return result.fetchall()

    def upsert_triplets(self, triplets: list[dict[str, int]]) -> None:
        with self.session() as session:
            session.execute(sqlite_upsert(Triplets).on_conflict_do_nothing(), triplets)
            session.commit()

    def upsert_references(self, references: list[dict[str, int]]) -> None:
        with self.session() as session:
            stmt = sqlite_upsert(References).values(references)
            session.execute(
                stmt.on_conflict_do_update(
                    stmt.table.primary_key, set_={"doi": stmt.excluded.doi}
                )
            )
            session.commit()

    def upsert_structures(self, structures: list[dict[str, int]]) -> None:
        with self.session() as session:
            for i in range(0, len(structures), 1000):
                stmt = sqlite_upsert(Structures).values(structures[i : i + 1000])
                session.execute(
                    stmt.on_conflict_do_update(
                        stmt.table.primary_key, set_={"smiles": stmt.excluded.smiles}
                    )
                )

                session.commit()

    def upsert_taxo_names(self, taxo_names: list[dict[str, int]]) -> None:
        with self.session() as session:
            for i in range(0, len(taxo_names), 1000):
                stmt = sqlite_upsert(TaxoNames).values(taxo_names[i : i + 1000])
                session.execute(
                    stmt.on_conflict_do_update(
                        stmt.table.primary_key, set_={"name": stmt.excluded.name}
                    )
                )
                session.commit()

    def upsert_rank_names(self, ranks_names: list[dict[str, int]]) -> None:
        with self.session() as session:
            for i in range(0, len(ranks_names), 1000):
                stmt = sqlite_upsert(TaxoRankNames).values(ranks_names[i : i + 1000])
                session.execute(
                    stmt.on_conflict_do_update(
                        stmt.table.primary_key, set_={"name": stmt.excluded.name}
                    )
                )
                session.commit()

    def upsert_taxo_ranks(self, taxo_ranks: list[dict[str, int]]) -> None:
        with self.session() as session:
            for i in range(0, len(taxo_ranks), 1000):
                stmt = sqlite_upsert(TaxoRanks).values(taxo_ranks[i : i + 1000])
                session.execute(
                    stmt.on_conflict_do_nothing() # TODO We probably want to generate a warning, or just wipe everything and rebuild
                )
                session.commit()

    def upsert_taxo_parenting(self, parenting: dict[int, dict[int, int]]) -> None:
        with self.session() as session:
            for item, parents in parenting.items():
                stmt = sqlite_upsert(TaxoParents).values(
                    [{"id": item, "parent_id": parent, "distance": distance} for parent, distance in parents.items()]
                )
                session.execute(
                    stmt.on_conflict_do_update(
                        stmt.table.primary_key, set_={"distance": stmt.excluded.distance}
                    )
                )
            session.commit()

    def get_generic_of_generic(self, out: Any, inp: Any, item: int) -> set[int]:
        with self.session() as session:
            result = session.query(out).filter(inp == item).distinct()
            return {row[0] for row in result}

    def get_generics_of_generics(self, out: Any, inp: Any, items: set[int]) -> set[int]:
        if len(items) > self.list_limit:
            raise Exception(
                f"Too many {inp}, the limit is {self.list_limit} and you have {len(items)}"
            )

        with self.session() as session:
            result = session.query(out).filter(inp.in_(items)).distinct()
            return {row[0] for row in result}

    def get_triplets_for(
        self,
        reference_ids: set[int] | None,
        structure_ids: set[int] | None,
        taxon_ids: set[int] | None,
    ) -> set[tuple[int, int, int]]:
        with self.session() as session:
            filters = []
            if reference_ids is not None:
                filters += [Triplets.reference_id.in_(reference_ids)]
            if structure_ids is not None:
                filters += [Triplets.structure_id.in_(structure_ids)]
            if taxon_ids is not None:
                filters += [Triplets.taxon_id.in_(taxon_ids)]
            result = session.query(
                Triplets.reference_id, Triplets.structure_id, Triplets.taxon_id
            )

            if len(filters) > 0:
                result = result.filter(*filters)

            return {(row[0], row[1], row[2]) for row in result}
