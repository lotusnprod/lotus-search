#!/usr/bin/env python3
import csv
import pickle


def run() -> None:
    database = {}

    with open("./data/database_chemo.pkl", "rb") as f:
        database.update(pickle.load(f))
    print("Finished integrating chemo")

    with open("./data/database_taxo.pkl", "rb") as f:
        database.update(pickle.load(f))
    print("Finished integrating taxo")

    # TODO references

    with open("./data/database.pkl", "wb") as f:
        pickle.dump(database, f)
    print("Finished dumping")


if __name__ == "__main__":
    run()
