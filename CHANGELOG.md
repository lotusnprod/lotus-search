# Changelog

## Unreleased

### 2023-12-27
    
    - More refactoring of SQLite
    - Move DOIs and SMILESes to the DB

### 2023-12-26

    - Dependencies update
    - Switch to SQLAlchemy
    - Move the schema of the db to version 2, delete your index.db
    - Triplet endpoints
    - Tests for API

### 2023-12-25

    - Simplify code
    - Switch to SQLite (WIP)

### 2023-12-17

    - Refactor code
    - Add delay when retrying queries
    - Remove dead code and projects
    - Added tests for the update script and fixtures
    - Fixed a bug where CSVs for taxonomy weren't loaded

### 2023-12-16

    - Fix tests
    - Refactor the download of CSVs
    - Parallelize the download of CSVs
    - Add parameters to the update script
        you can run with `--only xxx`, `--stop xxx` or `--skip xxx` to:
        - only execute xxx
        - stop before xxx
        - skip xxx 
        Where xxx can be a task or a group name
        Task/Groups descriptions can be listed with `--list`
    - Fix a potential abuse if taxon contains a comma in `generate_database_taxo.py`

## Release

### 0.1.2

