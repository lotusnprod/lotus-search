from update import (download_query_as_csv, generate_database,
                    generate_database_chemo, generate_database_index,
                    generate_database_taxo)
from update.models import Group, Task

DownloadGroup = Group(name="downloads", parallel=True)
DatabaseGroup = Group(name="database", parallel=True)
MergingGroup = Group(name="merging", parallel=False)


TASKS = [
    Task(
        name="references",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/references.rq",
            "output_file": "references.csv",
        },
    ),
    Task(
        name="ranks_names",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/ranks_names.rq",
            "output_file": "ranks_names.csv",
        },
    ),
    Task(
        name="structures",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/structures.rq",
            "output_file": "structures.csv",
        },
    ),
    Task(
        name="taxa",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={"query_file": "update/queries/taxa.rq", "output_file": "taxa.csv"},
    ),
    Task(
        name="taxa_names",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_names.rq",
            "output_file": "taxa_names.csv",
        },
    ),
    Task(
        name="taxa_parents",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_parents.rq",
            "output_file": "taxa_parents.csv",
        },
    ),
    Task(
        name="taxa_ranks",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_ranks.rq",
            "output_file": "taxa_ranks.csv",
        },
    ),
    Task(
        name="triplets",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/triplets.rq",
            "output_file": "triplets.csv",
        },
    ),
    Task(
        name="generate_database_chemo",
        f=generate_database_chemo.run,
        group=DatabaseGroup,
    ),
    Task(
        name="generate_database_taxo", f=generate_database_taxo.run, group=DatabaseGroup
    ),
    Task(
        name="generate_database_index",
        f=generate_database_index.run,
        group=MergingGroup,
    ),
    Task(name="generate_database", f=generate_database.run, group=MergingGroup),
]
MAX_WORKERS = 8
