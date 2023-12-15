#!/usr/bin/env python3

from update.common import run_query_to_csv, QLEVER_URL


def run() -> None:
    run_query_to_csv(
        query_file="update/queries/taxa_all.rq",
        output_file="taxa_all.csv",
        url = QLEVER_URL
        )
    

if __name__ == "__main__":
    run()
