import sys
import os

import polars as pl

CPG_ISLAND_CHR = "chr1"
CPG_ISLAND_POS = (207_321_199, 207_322_019)

CPG_ISLAND_EXPR = (pl.col("chr") == CPG_ISLAND_CHR) & (
    pl.col("pos").is_between(
        CPG_ISLAND_POS[0],
        CPG_ISLAND_POS[1],
    )
)


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

    print("region", "filename", "N", "X", sep="\t")
    for filename in os.listdir(INPUT_PATH):
        path = f"{INPUT_PATH}/{filename}"
        df = pl.read_csv(path, separator="\t")
        for region, expr in [("global", True), ("cd55_cgi", CPG_ISLAND_EXPR)]:
            filtered = df.filter(expr)
            print(
                region,
                filename,
                filtered["N"].sum(),
                filtered["X"].sum(),
                sep="\t",
            )
