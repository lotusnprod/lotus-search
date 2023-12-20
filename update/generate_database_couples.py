#!/usr/bin/env python3
import logging
import pickle
from csv import DictReader
from pathlib import Path

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def run(path: Path) -> None:
    t2s: dict[int, set[int]] = {}
    s2t: dict[int, set[int]] = {}
    ts2r: dict[tuple[int, int], set[str]] = {}

    with open(path / "couples.csv", "r") as f:
        reader = DictReader(f)
        for row in reader:
            s = int(row["structure"])
            t = int(row["taxon"])
            r = int(row["reference"])
            if t not in t2s:
                t2s[t] = set()
            if s not in s2t:
                s2t[s] = set()
            if (t, s) not in ts2r:
                ts2r[(t, s)] = set()
            t2s[t].add(s)
            s2t[s].add(t)
            ts2r[(t, s)].add(r)

    logging.info("Finished generating couples")

    database = {
        "t2s": t2s,
        "s2t": s2t,
        "ts2r": ts2r,
    }

    with open(path / "database_couples.pkl", "wb") as f:
        pickle.dump(database, f)
    logging.info("Finished dumping")


if __name__ == "__main__":
    run(Path("data"))
