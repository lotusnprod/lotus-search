# Changelog

## 0.1.2

### 2023-12-16
* Fix tests
* Refactor the download of CSVs
* Parallelize the download of CSVs
* Add parameters to the update script you can run with `--only xxx`, `--stop xxx` or `--skip xxx` to:
    * only execute xxx
    * stop before xxx
    * skip xxx
    * Where xxx can be a task or a group name Task/Groups descriptions can be listed with `--list`
* Fix a potential abuse if taxon contains a comma in `generate_database_taxo.py`
