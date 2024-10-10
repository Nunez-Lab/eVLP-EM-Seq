#!/usr/bin/env python3

import gzip
import sys


def usage():
    print(
        f"usage: {sys.argv[0]} <gzipped bismark bedGraph coverage (.cov) file> <output file for dss>"
    )


if len(sys.argv) == 1:
    usage()
    sys.exit(0)

if len(sys.argv) != 3:
    print("error: incorrect number of arguments")
    usage()
    sys.exit(1)

# See the following link for why we need this conversion:
# https://bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.html#3_Using_DSS_for_BS-seq_differential_methylation_analysis

# Should be gzip with columns: chr, start, end, methylation%, count methylated, count unmethylated
bismark_file = sys.argv[1]

# Will have columns: chr, start, total number of reads, number of reads showing methylation
output_file = sys.argv[2]

line_count = 0
with gzip.open(bismark_file, "r") as f:
    for _ in f:
        line_count += 1

with open(output_file, "w") as out_f:
    with gzip.open(bismark_file, "r") as in_f:
        out_f.write("chr\tpos\tN\tX\n")
        for i, line in enumerate(in_f):
            line = line.decode("utf-8")
            chromosome, start, _, _, methylated, unmethylated = (
                line.strip().split("\t")
            )
            methylated = int(methylated)
            unmethylated = int(unmethylated)
            total = methylated + unmethylated
            out_f.write(f"{chromosome}\t{start}\t{total}\t{methylated}\n")

            if i % 1000000 == 0:
                print(f"Completed {i} of {line_count} for {bismark_file}")
