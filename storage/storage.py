import sqlite3
from pathlib import Path
from typing import Any

from sqlalchemy import create_engine, insert, text
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
            self.create_tables()

    def session(self, autoflush=True):
        Session = sessionmaker(bind=self.engine, autoflush=autoflush)
        return Session()

    def query(self, query: str):
        with self.engine.connect() as connection:
            result = connection.execute(text(query))
        return result.fetchall()

    def upsert_triplets(self, triplets: list[dict[str, int]]) -> None:
        with self.session(autoflush=False) as session:
            for i in range(0, len(triplets), self.list_limit // 3):
                session.execute(
                    insert(Triplets), triplets[i : i + self.list_limit // 3]
                )
            session.commit()

    def upsert_references(self, references: list[dict[str, object]]) -> None:
        with self.session(autoflush=False) as session:
            for i in range(0, len(references), self.list_limit // 2):
                session.execute(
                    insert(References).values(),
                    references[i : i + self.list_limit // 2],
                )
            session.commit()

    def upsert_structures(self, structures: list[dict[str, object]]) -> None:
        with self.session(autoflush=False) as session:
            for i in range(0, len(structures), self.list_limit // 2):
                session.execute(
                    insert(Structures).values(structures[i : i + self.list_limit // 2])
                )
            session.commit()

    def upsert_taxo_names(self, taxo_names: list[dict[str, object]]) -> None:
        with self.session(autoflush=False) as session:
            for i in range(0, len(taxo_names), self.list_limit // 2):
                session.execute(
                    insert(TaxoNames), taxo_names[i : i + self.list_limit // 2]
                )
            session.commit()

    def upsert_rank_names(self, ranks_names: list[dict[str, object]]) -> None:
        with self.session(autoflush=False) as session:
            for i in range(0, len(ranks_names), self.list_limit // 2):
                session.execute(
                    insert(TaxoRankNames), ranks_names[i : i + self.list_limit // 2]
                )
            session.commit()

    def upsert_taxo_ranks(self, taxo_ranks: list[dict[str, int]]) -> None:
        with self.session(autoflush=False) as session:
            for i in range(0, len(taxo_ranks), self.list_limit // 2):
                session.execute(
                    insert(TaxoRanks), taxo_ranks[i : i + self.list_limit // 2]
                )
            session.commit()

    def upsert_taxo_parenting(self, parenting: list[tuple[int, int, int]]) -> None:
        with self.session(autoflush=False) as session:
            for i in range(0, len(parenting), self.list_limit // 2):
                session.execute(
                    insert(TaxoParents),
                    [
                        {"id": item[0], "parent_id": item[1], "distance": item[2]}
                        for item in parenting[i : i + self.list_limit // 2]
                    ],
                )

            session.commit()

    def get_generic_of_generic(self, out: Any, inp: Any, item: int) -> set[int]:
        with self.session() as session:
            result = session.query(out).filter(inp == item).distinct()
            return {row[0] for row in result}

    def get_generics_of_generics(self, out: Any, inp: Any, items: set[int]) -> set[int]:
        items_set = list(items)
        output = set()

        with self.session() as session:
            for i in range(0, len(items_set), self.list_limit):
                result = session.query(out).filter(
                    inp.in_(items_set[i : i + self.list_limit])
                )
                output |= {row[0] for row in result}

        return output

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

    def create_tables(self):
        Base.metadata.create_all(self.engine)
        with self.session() as session:
            session.add(SchemaVersion(version=self.SCHEMA_VERSION))
            session.commit()

    def drop_and_create_tables(self):
        Base.metadata.drop_all(self.engine)
        self.create_tables()
