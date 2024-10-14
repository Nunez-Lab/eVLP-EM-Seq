import sys

import polars as pl


def usage():
    print(
        "usage: {sys.argv[0]} <path to OT dss input> <path to OB dss input> <desired path to combined input>"
    )


if __name__ == "__main__":
    if len(sys.argv) == 1:
        usage()
        sys.exit(0)

    if len(sys.argv) != 4:
        print("error: incorrect number of arguments")
        usage()
        sys.exit(1)

    TOP_PATH = sys.argv[1]
    BOTTOM_PATH = sys.argv[2]
    OUTPUT_PATH = sys.argv[3]

    top = pl.scan_csv(TOP_PATH, separator="\t")
    bot = pl.scan_csv(BOTTOM_PATH, separator="\t")

    combined = (
        pl.concat([top, bot.with_columns(pos=pl.col("pos") - 1)])
        .group_by(["chr", "pos"])
        .sum()
        .collect()
    )

    combined.write_csv(OUTPUT_PATH, separator="\t")
