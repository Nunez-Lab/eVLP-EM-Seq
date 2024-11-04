# %% Import

import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import polars as pl
import scipy.stats as st

# %% Get command line arguments

# sys.argv = [sys.argv[0]]
# sys.argv.append("../data/dss-output/dss.csv")
# sys.argv.append(
#     "../rna-seq-correlation-data/deseq2_results_with_gene_names.csv"
# )
# sys.argv.append("../rna-seq-correlation-data/MANE.GRCh38.v1.4.summary.txt")
# sys.argv.append("../rna-seq-correlation-data/chr_names.csv")
# sys.argv.append("../output")

METHYL_PATH = sys.argv[1]
RNA_PATH = sys.argv[2]
TSS_PATH = sys.argv[3]
CHR_PATH = sys.argv[4]
OUTPUT_PATH = sys.argv[5]

# %% Load data

raw_methyl = pl.scan_csv(METHYL_PATH)
raw_rna = pl.scan_csv(RNA_PATH)
raw_mane = pl.scan_csv(TSS_PATH, separator="\t")
raw_chr = pl.scan_csv(CHR_PATH)

# %% Preprocess data

methyl_loci = (
    raw_methyl.filter(~pl.col("chr").str.contains("_"))
    .with_columns(
        chr=pl.when(pl.col("chr") == "chrX")
        .then(pl.lit("X"))
        .when(pl.col("chr") == "chrY")
        .then(pl.lit("Y"))
        .when(pl.col("chr") == "chrM")
        .then(pl.lit("MT"))
        .otherwise(pl.col("chr").str.slice(3)),
        pval_methyl=pl.col("pval"),
    )
    .join(raw_chr, on="chr")
    .select(
        "GRCh38_chr",
        "pos",
        "pval_methyl",
    )
)

rna_gene = (
    raw_rna.filter(pl.col("padj") != "NA")
    .with_columns(
        pval_rna=pl.col("pvalue").cast(pl.Float64),
        Ensembl_Gene=pl.col("Geneid"),
        symbol=pl.col("Gene.name"),
    )
    .select(
        "Ensembl_Gene",
        "symbol",
        "pval_rna",
    )
)

mane_gene = raw_mane.with_columns(
    promoter_start=pl.col("chr_start") - 500,
    promoter_end=pl.col("chr_start") + 500,
).select(
    "Ensembl_Gene",
    "GRCh38_chr",
    "promoter_start",
    "promoter_end",
)

# %% Annotate methylation data

methyl_gene = (
    methyl_loci.sort(by="pos")
    .join_asof(
        mane_gene.sort(by="promoter_start"),
        by="GRCh38_chr",
        left_on="pos",
        right_on="promoter_start",
        strategy="backward",
        coalesce=False,
    )
    .filter(~pl.col("promoter_end").is_null())
    .filter(pl.col("pos") <= pl.col("promoter_end"))
    .group_by("Ensembl_Gene")
    .agg(min_pval_methyl=pl.col("pval_methyl").min())
)

# %% Combine data and compute scores

data = (
    methyl_gene.join(
        rna_gene,
        on="Ensembl_Gene",
    )
    .with_columns(
        score_methyl=pl.col("min_pval_methyl").log10().neg(),
        score_rna=pl.col("pval_rna").log10().neg(),
    )
    .select("symbol", "score_methyl", "score_rna")
    .collect()
)

# %% Compute correlation

corr = st.spearmanr(data["score_methyl"], data["score_rna"])

# %% Make plot

# Source: https://personal.sron.nl/~pault/
BLUE = "#4477AA"
PURPLE = "#AA3377"

MAX_SCORE = 150
SCORE_TICKS_EVERY = 30
COLOR_CUTOFF = 20

assert (data["score_methyl"] < MAX_SCORE).all(), data["score_methyl"].max()
assert (data["score_rna"] < MAX_SCORE).all(), data["score_rna"].max()

colors = data.select(
    pl.when((data["score_methyl"] > 20) & (data["score_rna"] > 20))
    .then(pl.lit(PURPLE))
    .otherwise(pl.lit(BLUE))
)["literal"]

_, cd55_methyl, cd55_rna = data.filter(data["symbol"] == "CD55").row(0)

fig, ax = plt.subplots(1, 1, figsize=(5, 5))

ax.scatter(
    data["score_methyl"],
    data["score_rna"],
    c=colors,
    marker=".",
)

ax.annotate(
    "CD55",
    xy=(cd55_methyl, cd55_rna),
    xytext=(-3, 3),
    textcoords="offset points",
    ha="right",
    fontweight="bold",
    color=PURPLE,
)

ax.text(
    0.95,
    0.95,
    r"Spearman’s $\rho = $"
    + str(round(corr.statistic, 3))
    + r" ($p = $"
    + str(round(corr.pvalue, 3))
    + ")",
    transform=ax.transAxes,
    ha="right",
)

ax.set_xticks(np.arange(0, MAX_SCORE + 1, SCORE_TICKS_EVERY))
ax.set_yticks(np.arange(0, MAX_SCORE + 1, SCORE_TICKS_EVERY))
ax.set_xlim([0, MAX_SCORE])
ax.set_ylim([0, MAX_SCORE])

ax.set_xlabel(
    r"$\bf{Methylation}$ $\bf{significance}$ $\bf{score}$"
    + "\n"
    + r"$-\log_{10}($min. promoter $p$-value$)$"
)

ax.set_ylabel(
    r"$\bf{Expression}$ $\bf{significance}$ $\bf{score}$"
    + "\n"
    + r"$-\log_{10}(p$-value$)$"
)

ax.spines[["top", "right"]].set_visible(False)

ax.set_aspect("equal", adjustable="box")
fig.tight_layout()

fig.savefig(f"{OUTPUT_PATH}/corr.png", dpi=300)
fig.savefig(f"{OUTPUT_PATH}/corr.svg")
