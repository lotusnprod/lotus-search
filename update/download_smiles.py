#!/usr/bin/env python3

from update.common import run_query_to_csv


def run() -> None:
    run_query_to_csv(
        query_file="update/queries/structures.rq",
        output_file="smiles.csv"
        )
    

if __name__ == "__main__":
    run()
