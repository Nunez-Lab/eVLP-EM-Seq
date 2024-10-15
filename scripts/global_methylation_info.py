import sys
import os

import polars as pl


def usage():
    print("usage: {sys.argv[0]} <path to combined+filtered files>")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        usage()
        sys.exit(0)

    if len(sys.argv) != 2:
        print("error: incorrect number of arguments")
        usage()
        sys.exit(1)

    INPUT_PATH = sys.argv[1]

    print("filename", "N", "X", sep="\t")
    for filename in os.listdir(INPUT_PATH):
        path = f"{INPUT_PATH}/{filename}"
        df = pl.read_csv(path, separator="\t")
        print(filename, df["N"].sum(), df["X"].sum(), sep="\t")
