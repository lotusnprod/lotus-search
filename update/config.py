from update import (
    download_query_as_csv,
    generate_database_chemo,
    generate_database_index,
)
from update.models import Group, Task

DownloadGroup = Group(name="downloads", parallel=True)
DatabaseGroup = Group(name="database", parallel=True)
MergingGroup = Group(name="merging", parallel=False)


TASKS = [
    Task(
        name="lotus_enhanced",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/triplets_enhanced.rq",
            "url": "https://qlever.cs.uni-freiburg.de/api/wikidata",
            "output_file": "lotus_enhanced.csv",
        },
    ),
    Task(
        name="classes_cxsmiles",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/classes_cxsmiles.rq",
            "output_file": "classes_cxsmiles.csv",
        },
    ),
    Task(
        name="classes_parents",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/classes_parents.rq",
            "output_file": "classes_parents.csv",
        },
    ),
    Task(
        name="classes_smarts",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/classes_smarts.rq",
            "output_file": "classes_smarts.csv",
        },
    ),
    Task(
        name="classes_smiles_canonical",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/classes_smiles_canonical.rq",
            "output_file": "classes_smiles_canonical.csv",
        },
    ),
    Task(
        name="classes_smiles_isomeric",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/classes_smiles_isomeric.rq",
            "output_file": "classes_smiles_isomeric.csv",
        },
    ),
    Task(
        name="journals",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/journals.rq",
            "output_file": "journals.csv",
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
        name="references",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/references.rq",
            "output_file": "references.csv",
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
        name="structures_classes",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/structures_classes.rq",
            "output_file": "structures_classes.csv",
        },
    ),
    Task(
        name="structures_ids_chebi",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/structures_ids_chebi.rq",
            "output_file": "structures_ids_chebi.csv",
        },
    ),
    Task(
        name="structures_ids_chembl",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/structures_ids_chembl.rq",
            "output_file": "structures_ids_chembl.csv",
        },
    ),
    Task(
        name="structures_ids_chemspider",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/structures_ids_chemspider.rq",
            "output_file": "structures_ids_chemspider.csv",
        },
    ),
    Task(
        name="structures_ids_kegg",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/structures_ids_kegg.rq",
            "output_file": "structures_ids_kegg.csv",
        },
    ),
    Task(
        name="structures_ids_pubchem",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/structures_ids_pubchem.rq",
            "output_file": "structures_ids_pubchem.csv",
        },
    ),
    Task(
        name="structures_ids_zinc",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/structures_ids_zinc.rq",
            "output_file": "structures_ids_zinc.csv",
        },
    ),
    Task(
        name="taxa_ids_col",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_ids_col.rq",
            "output_file": "taxa_ids_col.csv",
        },
    ),
    Task(
        name="taxa_ids_eol",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_ids_eol.rq",
            "output_file": "taxa_ids_eol.csv",
        },
    ),
    Task(
        name="taxa_ids_gbif",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_ids_gbif.rq",
            "output_file": "taxa_ids_gbif.csv",
        },
    ),
    Task(
        name="taxa_ids_inat",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_ids_inat.rq",
            "output_file": "taxa_ids_inat.csv",
        },
    ),
    Task(
        name="taxa_ids_ipni",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_ids_ipni.rq",
            "output_file": "taxa_ids_ipni.csv",
        },
    ),
    Task(
        name="taxa_ids_irmng",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_ids_irmng.rq",
            "output_file": "taxa_ids_irmng.csv",
        },
    ),
    Task(
        name="taxa_ids_itis",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_ids_itis.rq",
            "output_file": "taxa_ids_itis.csv",
        },
    ),
    Task(
        name="taxa_ids_ncbi",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_ids_ncbi.rq",
            "output_file": "taxa_ids_ncbi.csv",
        },
    ),
    Task(
        name="taxa_ids_otl",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_ids_otl.rq",
            "output_file": "taxa_ids_otl.csv",
        },
    ),
    Task(
        name="taxa_ids_powo",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_ids_powo.rq",
            "output_file": "taxa_ids_powo.csv",
        },
    ),
    Task(
        name="taxa_ids_tropicos",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_ids_tropicos.rq",
            "output_file": "taxa_ids_tropicos.csv",
        },
    ),
    Task(
        name="taxa_ids_wfo",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_ids_wfo.rq",
            "output_file": "taxa_ids_wfo.csv",
        },
    ),
    Task(
        name="taxa_ids_worms",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_ids_worms.rq",
            "output_file": "taxa_ids_worms.csv",
        },
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
        name="taxa_names_com",
        f=download_query_as_csv.run,
        group=DownloadGroup,
        params={
            "query_file": "update/queries/taxa_names_com.rq",
            "output_file": "taxa_names_com.csv",
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
        name="generate_database_index",
        f=generate_database_index.run,
        group=MergingGroup,
    ),
]
MAX_WORKERS = 8
