import os
import sqlite3
import logging
import time
from pathlib import Path
from typing import Any

from sqlalchemy import create_engine, insert, text
from sqlalchemy.orm import sessionmaker

# Keep that this way so metadata gets all the tables
from storage.models import (
    Base,
    Journals,
    References,
    SchemaVersion,
    Structures,
    StructuresDescriptors,
    TaxoNames,
    TaxoParents,
    TaxoRankNames,
    TaxoRanks,
    Triplets,
)


class Storage:
    """Lightweight wrapper over SQLAlchemy engine/session for SQLite access."""

    SCHEMA_VERSION = 6  # Increment when schema changes (enforced in update pipeline)

    def __init__(self, path: Path):
        self.con = sqlite3.connect(path / "index.db")
        self.list_limit = self.con.getlimit(sqlite3.SQLITE_LIMIT_VARIABLE_NUMBER)
        self.con.close()
        self.db_path = str((path / "index.db").absolute())
        self.engine = create_engine(f"sqlite:///{self.db_path}")
        # Check if the schema table exists if not, call create
        new_db = self.query(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='schema_version'",
        )

        if len(new_db) == 0:
            self.create_tables()

    def session(self, autoflush: bool = True):
        """Create a new session."""
        Session = sessionmaker(bind=self.engine, autoflush=autoflush)
        return Session()

    def query(self, query: str):  # type: ignore[no-untyped-def]
        """Execute a query and return all results."""
        with self.engine.connect() as connection:
            result = connection.execute(text(query))
        return result.fetchall()

    def upsert_journals(self, journals: list[dict[str, object]]) -> None:
        """Upsert journals data."""
        with self.session(autoflush=False) as session:
            for i in range(0, len(journals), self.list_limit // 2):
                session.execute(
                    insert(Journals),
                    journals[i : i + self.list_limit // 2],
                )
            session.commit()

    def upsert_triplets(self, triplets: list[dict[str, int]]) -> None:
        """Upsert triplets data."""
        with self.session(autoflush=False) as session:
            for i in range(0, len(triplets), self.list_limit // 3):
                session.execute(
                    insert(Triplets),
                    triplets[i : i + self.list_limit // 3],
                )
            session.commit()

    def upsert_references(self, references: list[dict[str, object]]) -> None:
        """Upsert references data."""
        with self.session(autoflush=False) as session:
            for i in range(0, len(references), self.list_limit // 2):
                session.execute(
                    insert(References).values(),
                    references[i : i + self.list_limit // 2],
                )
            session.commit()

    def upsert_structures(self, structures: list[dict[str, object]]) -> None:
        """Upsert structures data."""
        with self.session(autoflush=False) as session:
            for i in range(0, len(structures), self.list_limit // 2):
                session.execute(
                    insert(Structures),
                    structures[i : i + self.list_limit // 2],
                )
            session.commit()

    def upsert_structures_descriptors(self, descriptors: dict[str, dict]) -> None:
        """Upsert structure descriptors safely respecting SQLite variable limits.

        Previous implementation could hit 'too many SQL variables' when the number
        of SMILES exceeded the backend parameter limit (seen as an OperationalError).
        We now batch both the lookup (IN clause) and the bulk insert.
        """
        if not descriptors:
            return
        # Allow override mainly for debugging / tuning
        default_lookup_batch = max(50, self.list_limit - 10)
        lookup_batch_size = int(
            os.getenv("LOTUS_DESCRIPTOR_LOOKUP_BATCH", default_lookup_batch)
        )
        # Each descriptor row has 3 bound parameters; keep comfortable margin
        default_insert_batch = max(50, (self.list_limit // 3) - 10)
        insert_batch_size = int(
            os.getenv("LOTUS_DESCRIPTOR_INSERT_BATCH", default_insert_batch)
        )

        with self.session(autoflush=False) as session:
            smiles_list = list(descriptors.keys())
            smiles_to_structure_id: dict[str, int] = {}
            for i in range(0, len(smiles_list), lookup_batch_size):
                batch = smiles_list[i : i + lookup_batch_size]
                if not batch:
                    continue
                for structure in (
                    session.query(Structures).filter(Structures.smiles.in_(batch)).all()
                ):
                    smiles_to_structure_id[structure.smiles] = structure.id

            descriptors_buffer: list[StructuresDescriptors] = []
            flush = session.bulk_save_objects

            for smiles, descriptor_data in descriptors.items():
                structure_id = smiles_to_structure_id.get(smiles)
                if structure_id is None:
                    continue
                for descriptor_name, descriptor_value in descriptor_data.items():
                    if descriptor_name == "smiles":
                        continue
                    descriptors_buffer.append(
                        StructuresDescriptors(
                            structure_id=structure_id,
                            descriptor_name=descriptor_name,
                            descriptor_value=descriptor_value,
                        ),
                    )
                    if len(descriptors_buffer) >= insert_batch_size:
                        flush(descriptors_buffer)
                        descriptors_buffer.clear()

            if descriptors_buffer:
                flush(descriptors_buffer)
            session.commit()

    def upsert_taxo_names(self, taxo_names: list[dict[str, object]]) -> None:
        """Upsert taxonomy names data."""
        with self.session(autoflush=False) as session:
            for i in range(0, len(taxo_names), self.list_limit // 2):
                session.execute(
                    insert(TaxoNames),
                    taxo_names[i : i + self.list_limit // 2],
                )
            session.commit()

    def upsert_rank_names(self, ranks_names: list[dict[str, object]]) -> None:
        """Upsert taxonomy rank names data."""
        with self.session(autoflush=False) as session:
            for i in range(0, len(ranks_names), self.list_limit // 2):
                session.execute(
                    insert(TaxoRankNames),
                    ranks_names[i : i + self.list_limit // 2],
                )
            session.commit()

    def upsert_taxo_ranks(self, taxo_ranks: list[dict[str, int]]) -> None:
        """Upsert taxonomy ranks data."""
        with self.session(autoflush=False) as session:
            for i in range(0, len(taxo_ranks), self.list_limit // 2):
                session.execute(
                    insert(TaxoRanks),
                    taxo_ranks[i : i + self.list_limit // 2],
                )
            session.commit()

    def upsert_taxo_parenting(self, parenting: list[tuple[int, int, int, int]]) -> None:
        """Upsert taxonomy parenting data (id, child_id, parent_id, distance).

        Enhancements (fast path):
        - Raw executemany batching to reduce ORM overhead.
        - Progress logging (rows & batches) controlled by env vars.
        - Optional post-insert index creation for faster downstream lookups.

        Environment variables:
        LOTUS_TAXO_PARENTING_FAST=1                 Use raw sqlite executemany path (default 1)
        LOTUS_TAXO_PARENTING_BATCH=NNN              Batch size (default adaptive ~50000)
        LOTUS_TAXO_PARENTING_PROGRESS=1             Enable progress logs (default 1)
        LOTUS_TAXO_PARENTING_PROGRESS_EVERY=N       Log every N batches (default 100)
        LOTUS_TAXO_PARENTING_INDEXES=1              Create helpful indexes after insert (default 1)
        LOTUS_TAXO_PARENTING_PRAGMAS=1              Apply speed PRAGMAs during load (default 1)
        LOTUS_TAXO_PARENTING_DROP_INDEXES=1         Drop existing indexes before load (default 1)
        LOTUS_TAXO_PARENTING_STATS_EVERY=N          Compute and log rolling r/s (default 20)
        """
        if not parenting:
            return

        fast = os.getenv("LOTUS_TAXO_PARENTING_FAST", "1") == "1"
        create_indexes = os.getenv("LOTUS_TAXO_PARENTING_INDEXES", "1") == "1"
        show_progress = os.getenv("LOTUS_TAXO_PARENTING_PROGRESS", "1") == "1"
        apply_pragmas = os.getenv("LOTUS_TAXO_PARENTING_PRAGMAS", "1") == "1"
        drop_indexes = os.getenv("LOTUS_TAXO_PARENTING_DROP_INDEXES", "1") == "1"
        progress_every = max(
            1, int(os.getenv("LOTUS_TAXO_PARENTING_PROGRESS_EVERY", "100"))
        )
        stats_every = max(1, int(os.getenv("LOTUS_TAXO_PARENTING_STATS_EVERY", "20")))

        total = len(parenting)

        # Determine batch size. Old heuristic tied to SQLITE variable limit (unnecessary for executemany).
        # Use large batches for fewer round-trips while keeping memory modest.
        if os.getenv("LOTUS_TAXO_PARENTING_BATCH"):
            batch_size = max(
                1000, int(os.getenv("LOTUS_TAXO_PARENTING_BATCH", "50000"))
            )
        else:
            # Adaptive: aim for about 2k batches (cap 100k) but at least 10k.
            target_batches = 2000
            batch_size = min(100000, max(10000, total // target_batches or 50000))
        # Safety upper bound: don't exceed 200k rows per executemany call.
        batch_size = min(batch_size, 200000)

        if not fast:
            # Original ORM chunked insert fallback (with progress logging)
            with self.session(autoflush=False) as session:
                batches = 0
                for i in range(0, total, self.list_limit // 2):
                    chunk = parenting[i : i + self.list_limit // 2]
                    session.execute(
                        insert(TaxoParents),
                        [
                            {
                                "id": item[0],
                                "child_id": item[1],
                                "parent_id": item[2],
                                "distance": item[3],
                            }
                            for item in chunk
                        ],
                    )
                    batches += 1
                    if show_progress and batches % progress_every == 0:
                        pct = (min(i + len(chunk), total) / total) * 100
                        logging.info(
                            "TAXO_INSERT orm batches=%d rows=%d/%d (%.2f%%)",
                            batches,
                            min(i + len(chunk), total),
                            total,
                            pct,
                        )
                session.commit()
        else:
            start = time.perf_counter()
            conn = sqlite3.connect(self.db_path)
            cur = conn.cursor()

            if apply_pragmas:
                # Fast import recommended PRAGMAs (volatile – okay for ephemeral bulk load).
                for pragma in (
                    "PRAGMA synchronous=OFF",
                    "PRAGMA journal_mode=MEMORY",
                    "PRAGMA temp_store=MEMORY",
                    "PRAGMA locking_mode=EXCLUSIVE",
                    "PRAGMA cache_size=-80000",  # ~80MB page cache
                ):
                    try:
                        cur.execute(pragma)
                    except Exception:  # pragma: no cover - defensive
                        pass

            if drop_indexes:
                # Drop per-row maintenance burden; recreate later.
                for idx in (
                    "taxo_parent_id",
                    "taxo_parent_child_id",
                    "taxo_parent_parent_id",
                    "taxo_parent_distance_id",
                    # Combined index might exist from previous runs
                    "idx_taxo_parents_child",
                    "idx_taxo_parents_parent",
                    "idx_taxo_parents_child_parent",
                ):
                    try:
                        cur.execute(f"DROP INDEX IF EXISTS {idx}")
                    except Exception:  # pragma: no cover
                        pass
                conn.commit()

            insert_sql = "INSERT INTO taxo_parents (id, child_id, parent_id, distance) VALUES (?,?,?,?)"
            batches = 0
            inserted = 0
            last_stat_inserted = 0
            last_stat_time = start
            progress_start_time = start

            # Wrap in explicit transaction to ensure PRAGMAs like synchronous=OFF are effective
            cur.execute("BEGIN TRANSACTION")
            for i in range(0, total, batch_size):
                chunk = parenting[i : i + batch_size]
                cur.executemany(insert_sql, chunk)
                inserted += len(chunk)
                batches += 1

                if show_progress and batches % progress_every == 0:
                    pct = (inserted / total) * 100
                    elapsed = time.perf_counter() - progress_start_time
                    logging.info(
                        "TAXO_INSERT fast batches=%d rows=%d/%d (%.2f%%) elapsed=%.1fs",
                        batches,
                        inserted,
                        total,
                        pct,
                        elapsed,
                    )
                # Rolling stats (rows/sec over recent window)
                if show_progress and batches % stats_every == 0:
                    now = time.perf_counter()
                    window_rows = inserted - last_stat_inserted
                    window_time = now - last_stat_time
                    if window_time > 0:
                        rps = window_rows / window_time
                        avg_rps = inserted / (now - start)
                        logging.info(
                            "TAXO_INSERT rate recent=%.0f r/s avg=%.0f r/s batch=%d size=%d",
                            rps,
                            avg_rps,
                            batches,
                            batch_size,
                        )
                    last_stat_inserted = inserted
                    last_stat_time = now

            conn.commit()
            if show_progress:
                elapsed = time.perf_counter() - start
                logging.info(
                    "TAXO_INSERT fast COMPLETE rows=%d total=%.1fs batch=%d avg_rps=%.0f",
                    inserted,
                    elapsed,
                    batch_size,
                    inserted / elapsed if elapsed else 0,
                )
            if create_indexes:
                idx_start = time.perf_counter()
                # Create helpful indexes (IF NOT EXISTS keeps it idempotent)
                for sql in (
                    "CREATE INDEX IF NOT EXISTS idx_taxo_parents_child ON taxo_parents(child_id)",
                    "CREATE INDEX IF NOT EXISTS idx_taxo_parents_parent ON taxo_parents(parent_id)",
                    "CREATE INDEX IF NOT EXISTS idx_taxo_parents_child_parent ON taxo_parents(child_id,parent_id)",
                    "CREATE INDEX IF NOT EXISTS idx_taxo_parents_distance ON taxo_parents(distance)",
                ):
                    cur.execute(sql)
                conn.commit()
                if show_progress:
                    logging.info(
                        "TAXO_INSERT indexes created in %.2fs",
                        time.perf_counter() - idx_start,
                    )
            conn.close()
        # No return value; side effects committed

    def get_generic_of_generic(self, out: Any, inp: Any, item: int) -> set[int]:
        """Get generic items of a generic item."""
        with self.session() as session:
            result = session.query(out).filter(inp == item).distinct()
            return {row[0] for row in result}

    def get_generics_of_generics(self, out: Any, inp: Any, items: set[int]) -> set[int]:
        """Get generics of multiple generic items."""
        items_set = list(items)
        output: set[int] = set()

        with self.session() as session:
            for i in range(0, len(items_set), self.list_limit):
                result = session.query(out).filter(
                    inp.in_(items_set[i : i + self.list_limit]),
                )
                output |= {row[0] for row in result}

        return output

    def get_triplets_for(
        self,
        reference_ids: set[int] | None,
        structure_ids: set[int] | None,
        taxon_ids: set[int] | None,
    ) -> set[tuple[int, int, int]]:
        """Get triplets for given reference, structure, and taxon IDs."""
        with self.session() as session:
            filters = []
            if reference_ids is not None:
                filters += [Triplets.reference_id.in_(reference_ids)]
            if structure_ids is not None:
                filters += [Triplets.structure_id.in_(structure_ids)]
            if taxon_ids is not None:
                filters += [Triplets.taxon_id.in_(taxon_ids)]
            result = session.query(
                Triplets.reference_id,
                Triplets.structure_id,
                Triplets.taxon_id,
            )

            if len(filters) > 0:
                result = result.filter(*filters)

            return {(row[0], row[1], row[2]) for row in result}

    def create_tables(self):
        """Create database tables and initial schema version entry."""
        Base.metadata.create_all(self.engine)
        with self.session() as session:
            session.add(SchemaVersion(version=self.SCHEMA_VERSION))
            session.commit()

    def drop_and_create_tables(self):
        """Drop and recreate all database tables."""
        Base.metadata.drop_all(self.engine)
        self.create_tables()
