#!/usr/bin/env python3

from update.common import run_query_to_csv


def run() -> None:
    run_query_to_csv(
        query_file="update/queries/taxa.rq",
        output_file="taxa.csv"
        )
    

if __name__ == "__main__":
    run()
