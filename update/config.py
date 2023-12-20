from update import (download_query_as_csv, generate_database,
                    generate_database_biblio, generate_database_chemo,
                    generate_database_couples, generate_database_taxo)
from update.common import QLEVER_URL
from update.models import Group, Task

DownloadGroup = Group(name="downloads", parallel=True)
DatabaseGroup = Group(name="database", parallel=True)
MergingGroup = Group(name="merging", parallel=False)


TASKS = [
    Task(
        name="couples_referenced",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/couples_referenced.rq",
            "output_file": "couples.csv",
        },
    ),
    Task(
        name="dois",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/references.rq",
            "output_file": "dois.csv",
        },
    ),
    Task(
        name="smiles",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/structures.rq",
            "output_file": "smiles.csv",
        },
    ),
    Task(
        name="taxa",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={"query_file": "update/queries/taxa.rq", "output_file": "taxa.csv"},
    ),
    Task(
        name="taxa_all",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_all.rq",
            "output_file": "taxa_all.csv",
            "url": QLEVER_URL,
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
        name="generate_database_couples",
        f=generate_database_couples.run,
        group=DatabaseGroup,
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
        name="generate_database_biblio",
        f=generate_database_biblio.run,
        group=DatabaseGroup,
    ),
    Task(name="generate_database", f=generate_database.run, group=MergingGroup),
]
MAX_WORKERS = 7
