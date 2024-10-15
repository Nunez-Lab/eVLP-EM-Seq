# %% Import

import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import polars as pl

# %% Get command line arguments

GLOBAL_METHYLATION_INFO_PATH = sys.argv[1]
DSS_RESULTS_PATH = sys.argv[2]
OUTPUT = sys.argv[3]

# %% Load data

raw_methyl = pl.read_csv(GLOBAL_METHYLATION_INFO_PATH, separator="\t")
raw_data = pl.read_csv(DSS_RESULTS_PATH)

# %% Process data
# Important: Assumes mu1 is CD55off, mu2 is WT

data = (
    raw_data.with_columns(
        cd55off_mu=pl.col("mu1"),
        wt_mu=pl.col("mu2"),
        effect_size=pl.col("diff"),
        score=pl.when(pl.col("mu1") > pl.col("mu2"))
        .then(pl.col("pval").log10().neg())
        .otherwise(pl.col("pval").log10()),
        chr=pl.when(pl.col("chr") == "phage_lambda")
        .then(pl.col("chr"))
        .when(pl.col("chr") == "plasmid_puc19c")
        .then(pl.col("chr"))
        .when(pl.col("chr").str.contains("_"))
        .then(pl.lit("other"))
        .otherwise(pl.col("chr")),
    )
    .with_columns(
        chr_order=pl.when(pl.col("chr") == "chrX")
        .then(23)
        .when(pl.col("chr") == "chrY")
        .then(24)
        .when(pl.col("chr") == "chrM")
        .then(25)
        .when(pl.col("chr") == "other")
        .then(26)
        .when(pl.col("chr") == "phage_lambda")
        .then(27)
        .when(pl.col("chr") == "plasmid_puc19c")
        .then(28)
        .otherwise(pl.col("chr").str.slice(3).str.to_integer(strict=False)),
        nice_chr=pl.when(pl.col("chr") == "phage_lambda")
        .then(pl.lit("Unmethylated control\n(lambda phage DNA)"))
        .when(pl.col("chr") == "plasmid_puc19c")
        .then(pl.lit("Methylated control\n(pUC19 plasmid DNA)"))
        .otherwise(pl.col("chr")),
        short_chr=pl.when(pl.col("chr") == "other")
        .then(pl.lit("?"))
        .when(pl.col("chr") == "phage_lambda")
        .then(pl.lit("Lambda"))
        .when(pl.col("chr") == "plasmid_puc19c")
        .then(pl.lit("pUC19"))
        .otherwise(pl.col("chr").str.slice(3)),
    )
    .select(
        "chr",
        "nice_chr",
        "short_chr",
        "chr_order",
        "pos",
        "cd55off_mu",
        "wt_mu",
        "score",
        "effect_size",
    )
    .sort(by=["chr_order", "pos"])
    .with_row_index()
)

# %% Highest and lowest summary

summary_cols = ["chr", "pos", "cd55off_mu", "wt_mu", "score", "effect_size"]
data[summary_cols].top_k(100, by="score").write_csv(
    f"{OUTPUT}/cd55off_top_hits.tsv", separator="\t"
)
data[summary_cols].top_k(100, by="score", reverse=True).write_csv(
    f"{OUTPUT}/wt_top_hits.tsv", separator="\t"
)

# %% Manhattan plot helper

BLUE = "#4477AA"
RED = "#AA3377"
GREEN = "#228833"


def manhattan(
    df,
    *,
    by,
    feature,
    two_sided,
    highlight,
    index="index",
    main_color=BLUE,
    secondary_color=RED,
    highlight_color=GREEN,
    yticks=None,
    yticklabels=None,
    ylabel=None,
    yaxis_formatter=None,
    xtick_rotation=90,
    xlabel=None,
    figsize=(10, 3),
    use_xticks=True,
):
    if yticks is not None:
        assert df[feature].min() > min(yticks), df[feature].min()
        assert df[feature].max() < max(yticks), df[feature].max()

    if yticklabels is not None:
        assert yticks is not None

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    xticks = []
    xticklabels = []

    for (c,), g in df.group_by(by):
        lo = g[index].min()
        hi = g[index].max()
        mid = (lo + hi) / 2

        if use_xticks:
            xticks.append(mid)
            xticklabels.append(c)

        if two_sided is not None:
            g = g.with_columns(
                color=pl.when(pl.col(feature) > 0)
                .then(pl.lit(secondary_color))
                .otherwise(pl.lit(main_color))
            )
        else:
            g = g.with_columns(color=pl.lit(main_color))

        if highlight is not None:
            g = g.with_columns(
                color=pl.when(highlight)
                .then(pl.lit(highlight_color))
                .otherwise(pl.col("color"))
            )

        ax.scatter(
            g[index],
            g[feature],
            marker=",",
            s=1,
            c=g["color"],
        )

        if use_xticks:
            ax.axvline(
                lo,
                color="lightgray",
                linewidth=0.5,
            )

    left = df[index].min()
    right = df[index].max()

    if use_xticks:
        ax.axvline(right, color="lightgray", linewidth=0.5)

    ax.set_xlim(left, right)

    if xlabel:
        ax.set_xlabel(xlabel)

    ax.set_xticks(
        xticks,
        labels=xticklabels,
        rotation=xtick_rotation,
    )

    ax.set_ylabel(ylabel)

    if two_sided is not None:
        xpadding = 0.005
        ypadding = 0.02
        ax.text(
            1 - xpadding,
            1 - ypadding,
            f"Higher in {two_sided[0]}",
            color=secondary_color,
            fontweight="bold",
            ha="right",
            va="top",
            transform=ax.transAxes,
        )
        ax.text(
            1 - xpadding,
            0 + ypadding,
            f"Higher in {two_sided[1]}",
            color=main_color,
            fontweight="bold",
            ha="right",
            va="bottom",
            transform=ax.transAxes,
        )

    if yticks is not None:
        ax.set_yticks(
            yticks,
            labels=yticklabels if yticklabels is not None else yticks,
        )
    elif two_sided is not None:
        ylim_lo, ylim_hi = ax.get_ylim()
        bound = max(abs(ylim_lo), abs(ylim_hi))
        ax.set_ylim(-bound, bound)

    if yaxis_formatter is not None:
        ax.yaxis.set_major_formatter(yaxis_formatter)

    fig.tight_layout()

    return fig, ax


# %% Plot controls

for condition in ["wt", "cd55off"]:
    feature = f"{condition}_mu"
    fig, ax = manhattan(
        data.filter(
            pl.col("chr_order") > 26,
        ),
        by="nice_chr",
        feature=feature,
        two_sided=None,
        highlight=pl.col(feature) > 0.5,
        highlight_color=RED,
        yticks=np.arange(0, 1.1, 0.2),
        ylabel=r"$\bf{\%}$ $\bf{reads}$ $\bf{methylated}$",
        yaxis_formatter=mtick.PercentFormatter(1.0),
        xtick_rotation=0,
    )

    fig.savefig(f"{OUTPUT}/controls_{condition}.png", dpi=300)
    plt.close(fig)

# %% Define important genomic loci

CD55_GUIDE_CHR = "chr1"
CD55_GUIDE_POS = (207_321_714, 207_321_735)
CD55_GUIDE_EXPR = (pl.col("chr") == CD55_GUIDE_CHR) & (
    pl.col("pos").is_between(
        CD55_GUIDE_POS[0],
        CD55_GUIDE_POS[1],
    )
)
CD55_GUIDE_START_INDEX = data.filter(CD55_GUIDE_EXPR)["index"].min()
CD55_GUIDE_END_INDEX = data.filter(CD55_GUIDE_EXPR)["index"].max()

# Source: https://www.uniprot.org/uniprotkb/P08174/genomic-coordinates

CD55_CHR = "chr1"
CD55_POS = (207_321_766, 207_359_607)

CD55_EXPR = (pl.col("chr") == CD55_CHR) & (
    pl.col("pos").is_between(
        CD55_POS[0],
        CD55_POS[1],
    )
)
CD55_SURROUNDING_EXPR = (pl.col("chr") == CD55_CHR) & (
    pl.col("pos").is_between(
        CD55_POS[0] - 10_000,
        CD55_POS[1] + 10_000,
    )
)
CD55_TSS_INDEX = data.filter(CD55_EXPR)["index"].min()

# Source: https://www.uniprot.org/uniprotkb/O43296/genomic-coordinates

ZN264_CHR = "chr19"
ZN264_POS = (57_191_914, 57_212_978)

ZN264_EXPR = (pl.col("chr") == ZN264_CHR) & (
    pl.col("pos").is_between(
        ZN264_POS[0],
        ZN264_POS[1],
    )
)
ZN264_SURROUNDING_EXPR = (pl.col("chr") == ZN264_CHR) & (
    pl.col("pos").is_between(
        ZN264_POS[0] - 10_000,
        ZN264_POS[1] + 10_000,
    )
)
ZN264_TSS_INDEX = data.filter(ZN264_EXPR)["index"].min()

# %% Plot differential methylation

score_yticks = np.arange(-160, 161, 40)
effect_size_yticks = np.arange(-1, 1.1, 0.25)
sample = False

for filter_expr, use_xticks, base_filename in [
    (pl.col("chr_order") < 25, True, "differential_methylation_all"),
    (CD55_SURROUNDING_EXPR, False, "differential_methylation_cd55"),
    (CD55_SURROUNDING_EXPR, False, "differential_methylation_cd55_sgRNA"),
    (ZN264_SURROUNDING_EXPR, False, "differential_methylation_zn264"),
]:
    for (
        feature,
        ylabel,
        yticks,
        yticklabels,
    ) in [
        (
            "score",
            r"$\bf{Significance}$ $\bf{score}$"
            + "\n"
            + r"$-\log_{10}(p$-value$)$",
            score_yticks,
            abs(score_yticks),
        ),
        (
            "effect_size",
            r"$\bf{Effect}$ $\bf{size}$",
            effect_size_yticks,
            [round(yt, 2) for yt in effect_size_yticks],
        ),
    ]:
        df = data.sample(10_000) if sample else data

        fig, ax = manhattan(
            df.filter(filter_expr),
            by="chr",
            feature=feature,
            two_sided=("CRISPRoff-eVLP treated", "untreated"),
            highlight=None,
            yticks=yticks,
            yticklabels=yticklabels,
            ylabel=ylabel,
            use_xticks=use_xticks,
        )

        if "zn264" in base_filename:
            text = "ZN264"
            xy = (ZN264_TSS_INDEX, 1.01)
            xytext = (ZN264_TSS_INDEX, 1.12)
        else:
            text = "CD55"
            xy = (CD55_TSS_INDEX, 1.01)
            xytext = (CD55_TSS_INDEX, 1.12)

        ax.annotate(
            text + (" [SUBSAMPLED DATA]" if sample else ""),
            xy=xy,
            xytext=xytext,
            ha="center",
            va="bottom",
            fontweight="bold",
            xycoords=ax.get_xaxis_transform(),
            textcoords=ax.get_xaxis_transform(),
            arrowprops=dict(
                facecolor="black",
                headwidth=8,
                headlength=8,
            ),
        )

        if "cd55_sgRNA" in base_filename:
            ax.hlines(
                [max(yticks)],
                [CD55_GUIDE_START_INDEX],
                [CD55_GUIDE_END_INDEX],
                color="red",
                linewidth=3,
                label="sgRNA",
            )
            ax.text(
                CD55_GUIDE_START_INDEX,
                max(yticks) - 20,
                "sgRNA",
                color="red",
            )
            ax.set_yticks(yticks)

        suffix = "-SUBSAMPLED" if sample else ""

        fig.tight_layout()
        fig.savefig(f"{OUTPUT}/{base_filename}-{feature}{suffix}.png", dpi=300)
        plt.close(fig)

# %% Genomic position plotter

fig, ax = manhattan(
    data.filter(CD55_SURROUNDING_EXPR),
    by="chr",
    feature="score",
    index="pos",
    two_sided=("CRISPRoff-eVLP treated", "untreated"),
    highlight=None,
    yticks=score_yticks,
    yticklabels=abs(score_yticks),
    ylabel=r"$\bf{Significance}$ $\bf{score}$",
    use_xticks=False,
    figsize=(15, 3),
)

coords = np.arange(
    CD55_POS[0] - 10_000,
    CD55_POS[1] + 10_000,
    1_000,
)
ax.set_xticks(
    coords,
    labels=coords - CD55_POS[0],
    rotation=90,
)
ax.set_xlabel("Distance from CD55 TSS")
fig.tight_layout()
fig.savefig(f"{OUTPUT}/test1.pdf")

CLOSE_CD55_SURROUNDING_EXPR = (pl.col("chr") == CD55_CHR) & (
    pl.col("pos").is_between(
        CD55_POS[0] - 2_000,
        CD55_POS[0] + 2_000,
    )
)

fig, ax = manhattan(
    data.filter(CLOSE_CD55_SURROUNDING_EXPR),
    by="chr",
    feature="score",
    index="pos",
    two_sided=("CRISPRoff-eVLP treated", "untreated"),
    highlight=None,
    yticks=score_yticks,
    yticklabels=abs(score_yticks),
    ylabel=r"$\bf{Significance}$ $\bf{score}$",
    use_xticks=False,
    figsize=(15, 3),
)

coords = np.arange(
    CD55_POS[0] - 2_000,
    CD55_POS[0] + 2_000,
    100,
)
ax.set_xticks(
    coords,
    labels=coords - CD55_POS[0],
    rotation=90,
)
ax.set_xlabel("Distance from CD55 TSS")

CPG_ISLAND_CHR = "chr1"
CPG_ISLAND_POS = (207_321_199, 207_322_019)

ax.hlines(
    [max(score_yticks) + 3],
    [CPG_ISLAND_POS[0]],
    [CPG_ISLAND_POS[1]],
    color="red",
    linewidth=3,
)
ax.text(
    CPG_ISLAND_POS[0],
    max(score_yticks) - 12,
    "Annotated CGI",
    color="red",
)

fig.tight_layout()
fig.savefig(f"{OUTPUT}/test1.pdf")
plt.close(fig)

# %% Global methylation information

CONDITIONS = ["Untreated", "Treated\n(CRISPRoff\neVLP)"]
CONDITION_COLORS = [BLUE, RED]

methyl = (
    raw_methyl.with_columns(
        condition=pl.when(pl.col("filename").str.contains("WT_"))
        .then(0)
        .when(pl.col("filename").str.contains("CD55off_"))
        .then(1)
        .otherwise(9999),
        percent=pl.col("X") / pl.col("N"),
    )
    .group_by("condition")
    .map_groups(lambda df: df.with_row_index(name="replicate"))
)

methyl_agg = methyl[["condition", "percent"]].group_by("condition").mean()

fig, ax = plt.subplots(1, 1, figsize=(2.5, 3))

b = ax.bar(
    methyl_agg["condition"],
    methyl_agg["percent"],
    color=[CONDITION_COLORS[i] for i in methyl_agg["condition"]],
)

ax.bar_label(
    b,
    fmt=lambda x: str(round(x * 100)) + "%",
    padding=7,
    color="gray",
)

ax.scatter(
    methyl["condition"] + (methyl["replicate"] - 0.5) / 3,
    methyl["percent"],
    color="black",
)

ax.set_xlim(-0.5, 1.5)
ax.set_xticks([0, 1], labels=CONDITIONS)

ax.set_ylim(0, 1)
ax.set_ylabel("Average CpG methylation")
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

ax.spines[["top", "right"]].set_visible(False)

fig.tight_layout()
fig.savefig(f"{OUTPUT}/global_methylation.pdf")
plt.close(fig)
