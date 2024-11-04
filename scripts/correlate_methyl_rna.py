# %% Import

import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import polars as pl

# %% Get command line arguments

sys.argv = [sys.argv[0]]
sys.argv.append("../data/dss-output/dss.csv")
sys.argv.append(
    "../rna-seq-correlation-data/deseq2_results_with_gene_names.csv"
)
sys.argv.append("../rna-seq-correlation-data/MANE.GRCh38.v1.4.summary.txt")
sys.argv.append("../rna-seq-correlation-data/chr_names.csv")
sys.argv.append("../output")

METHYL_PATH = sys.argv[1]
RNA_PATH = sys.argv[2]
TSS_PATH = sys.argv[3]
CHR_PATH = sys.argv[4]
OUTPUT_PATH = sys.argv[5]

# %% Load data

raw_methyl = pl.read_csv(METHYL_PATH)
raw_rna = pl.read_csv(RNA_PATH)
raw_mane = pl.read_csv(TSS_PATH, separator="\t")
raw_chr = pl.read_csv(CHR_PATH)

# %% Preprocess data

methyl = (
    raw_methyl.filter(~pl.col("chr").str.contains("_"))
    .with_columns(
        chr=pl.when(pl.col("chr") == "chrX")
        .then(pl.lit("X"))
        .when(pl.col("chr") == "chrY")
        .then(pl.lit("Y"))
        .when(pl.col("chr") == "chrM")
        .then(pl.lit("MT"))
        .otherwise(pl.col("chr").str.slice(3)),
    )
    .join(raw_chr, on="chr")
    .select(
        "GRCh38_chr",
        "pos",
        "fdr",
    )
)

rna = (
    raw_rna.filter(pl.col("padj") != "NA")
    .with_columns(
        fdr=pl.col("padj").cast(pl.Float64),
        Ensembl_Gene=pl.col("Geneid"),
    )
    .select(
        "Ensembl_Gene",
        "fdr",
    )
)

mane = raw_mane.select(
    "Ensembl_Gene", "symbol", "GRCh38_chr", "chr_start", "chr_end"
)

# %% Combine data

methyl_gene = methyl.join_where(
    mane,
    pl.col("GRCh38_chr") == pl.col("GRCh38_chr_right"),
    pl.col("chr_start") <= pl.col("pos"),
    pl.col("pos") <= pl.col("chr_end"),
)
