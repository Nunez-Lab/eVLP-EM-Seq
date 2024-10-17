# %% Import

import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import polars as pl

# %% Get command line arguments

sys.argv = [sys.argv[0]]
sys.argv.append("../data/average-methylation-info/info.tsv")
sys.argv.append("../data/dss-output/dss.csv")
sys.argv.append("../output")

AVERAGE_METHYLATION_INFO_PATH = sys.argv[1]
DSS_RESULTS_PATH = sys.argv[2]
OUTPUT = sys.argv[3]

# %% Load data

raw_methyl = pl.read_csv(AVERAGE_METHYLATION_INFO_PATH, separator="\t")
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

# Source: https://personal.sron.nl/~pault/
BLUE = "#4477AA"
PURPLE = "#AA3377"
GREEN = "#228833"
RED = "#EE6677"


def manhattan(
    df,
    *,
    by,
    feature,
    two_sided,
    highlight,
    index="index",
    main_color=BLUE,
    secondary_color=PURPLE,
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

    ax.spines[["top", "right"]].set_visible(False)

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
        highlight_color=PURPLE,
        yticks=np.arange(0, 1.1, 0.2),
        ylabel=r"$\bf{\%}$ $\bf{reads}$ $\bf{methylated}$",
        yaxis_formatter=mtick.PercentFormatter(1.0),
        xtick_rotation=0,
    )

    fig.savefig(f"{OUTPUT}/controls_{condition}.svg")
    plt.close(fig)

# %% Define important genomic loci

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

# %% Differential methylation plot parameters

SCORE_YTICKS = np.arange(-160, 161, 40)
EFFECT_SIZE_YTICKS = np.arange(-1, 1.1, 0.25)
SAMPLE = False

# %% Plot differential methylation

for filter_expr, use_xticks, base_filename in [
    (pl.col("chr_order") < 25, True, "differential_methylation_all"),
    (CD55_SURROUNDING_EXPR, False, "differential_methylation_cd55"),
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
            SCORE_YTICKS,
            abs(SCORE_YTICKS),
        ),
        (
            "effect_size",
            r"$\bf{Effect}$ $\bf{size}$",
            EFFECT_SIZE_YTICKS,
            [round(yt, 2) for yt in EFFECT_SIZE_YTICKS],
        ),
    ]:
        df = data.sample(10_000) if SAMPLE else data

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
            text + (" [SUBSAMPLED DATA]" if SAMPLE else ""),
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

        suffix = "-SUBSAMPLED" if SAMPLE else ""

        fig.tight_layout()
        fig.savefig(f"{OUTPUT}/{base_filename}-{feature}{suffix}.png", dpi=300)
        plt.close(fig)

# %% Important loci for genomic (zoomed) position plot

CD55_TSS_ZOOM_OFFSET = 2_000

CD55_ZOOMED_EXPR = (pl.col("chr") == CD55_CHR) & (
    pl.col("pos").is_between(
        CD55_POS[0] - CD55_TSS_ZOOM_OFFSET,
        CD55_POS[0] + CD55_TSS_ZOOM_OFFSET,
    )
)

CD55_ZOOM_COORDS = np.arange(
    CD55_POS[0] - CD55_TSS_ZOOM_OFFSET,
    CD55_POS[0] + CD55_TSS_ZOOM_OFFSET + 1,
    200,
)

# Source: UCSC Genome Browser
CPG_ISLAND_CHR = "chr1"
CPG_ISLAND_POS = (207_321_199, 207_322_019)

CD55_GUIDE_CHR = "chr1"
CD55_GUIDE_POS = (207_321_714, 207_321_735)

# %% Make zoomed plot v1

fig, ax = manhattan(
    data.filter(CD55_ZOOMED_EXPR),
    by="chr",
    feature="score",
    index="pos",
    two_sided=("CRISPRoff-eVLP treated", "untreated"),
    highlight=None,
    yticks=SCORE_YTICKS,
    yticklabels=abs(SCORE_YTICKS),
    ylabel=r"$\bf{Significance}$ $\bf{score}$"
    + "\n"
    + r"$-\log_{10}(p$-value$)$",
    use_xticks=False,
    figsize=(6, 3),
)

ax.set_xticks(
    CD55_ZOOM_COORDS,
    labels=CD55_ZOOM_COORDS - CD55_POS[0],
    rotation=90,
)

ax.set_xlabel("Distance from CD55 TSS", fontweight="bold")

ax.set_yticks(SCORE_YTICKS)
ax.set_ylim(SCORE_YTICKS[0], SCORE_YTICKS[-1])

ax.hlines(
    [-30],
    [CPG_ISLAND_POS[0]],
    [CPG_ISLAND_POS[1]],
    color=GREEN,
    linewidth=3,
)

ax.text(
    CPG_ISLAND_POS[0],
    -40,
    "Annotated CGI",
    color=GREEN,
    ha="left",
    va="top",
)

ax.hlines(
    [-85],
    [CD55_GUIDE_POS[0]],
    [CD55_GUIDE_POS[1]],
    color=RED,
    linewidth=3,
)

ax.text(
    CD55_GUIDE_POS[0],
    -95,
    "sgRNA",
    color=RED,
    ha="left",
    va="top",
)

fig.tight_layout()
fig.savefig(f"{OUTPUT}/cd55_zoom_v1.svg")
plt.close(fig)

# %% Make zoomed plot v2

df = data.filter(CD55_ZOOMED_EXPR)

fig, ax = plt.subplots(1, 1, figsize=(7, 3))

ax.scatter(
    df["pos"],
    df["cd55off_mu"],
    c=PURPLE,
    marker="o",
    s=7,
    label="CRISPRoff-eVLP treated",
)
ax.scatter(
    df["pos"],
    df["wt_mu"],
    c=BLUE,
    marker="^",
    s=7,
    label="Untreated",
)

ax.set_xticks(
    CD55_ZOOM_COORDS,
    labels=CD55_ZOOM_COORDS - CD55_POS[0],
    rotation=90,
)
ax.set_xlabel("Distance from CD55 TSS", fontweight="bold")

ax.set_yticks(np.arange(0, 1.01, 0.2))
ax.set_ylim(0, 1)
ax.set_ylabel("% Methylated reads", fontweight="bold")
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

ax.hlines(
    [0.95],
    [CPG_ISLAND_POS[0]],
    [CPG_ISLAND_POS[1]],
    color=GREEN,
    linewidth=3,
)

ax.text(
    CPG_ISLAND_POS[0],
    0.92,
    "Annotated CGI",
    color=GREEN,
    ha="left",
    va="top",
)

ax.hlines(
    [0.80],
    [CD55_GUIDE_POS[0]],
    [CD55_GUIDE_POS[1]],
    color=RED,
    linewidth=3,
)

ax.text(
    CD55_GUIDE_POS[0],
    0.77,
    "sgRNA",
    color=RED,
    ha="left",
    va="top",
)

leg = ax.legend(
    loc="lower right",
    handletextpad=0,
    labelcolor="linecolor",
    prop=dict(weight="bold"),
    borderpad=0.5,
)

leg.get_frame().set_facecolor("0.95")
leg.get_frame().set_linewidth(0.0)

ax.spines[["top", "right"]].set_visible(False)

fig.tight_layout()
fig.savefig(f"{OUTPUT}/cd55_zoom_v2.svg", dpi=300)
plt.close(fig)

# %% Average methylation information

CONDITIONS = ["Untreated", "Treated\n(CRISPRoff\neVLP)"]
CONDITION_COLORS = [BLUE, PURPLE]

methyl = (
    raw_methyl.with_columns(
        condition=pl.when(pl.col("filename").str.contains("WT_"))
        .then(0)
        .when(pl.col("filename").str.contains("CD55off_"))
        .then(1)
        .otherwise(9999),
        region=pl.when(pl.col("region") == "global")
        .then(pl.lit("Global Methylation"))
        .when(pl.col("region") == "cd55_cgi")
        .then(pl.lit("CD55 CGI Methylation"))
        .otherwise(pl.lit("UNKNOWN REGION")),
        percent=pl.col("X") / pl.col("N"),
    )
    .group_by("region", "condition")
    .map_groups(lambda df: df.with_row_index(name="replicate"))
)

for (region,), df in methyl.group_by("region"):
    fig, ax = plt.subplots(1, 1, figsize=(2.5, 3))

    agg = df[["condition", "percent"]].group_by("condition").mean()
    b = ax.bar(
        agg["condition"],
        agg["percent"],
        color=[CONDITION_COLORS[i] for i in agg["condition"]],
    )

    ax.bar_label(
        b,
        fmt=lambda x: str(round(x * 100)) + "%",
        padding=7,
        color="gray",
    )

    ax.scatter(
        df["condition"] + (df["replicate"] - 0.5) / 3,
        df["percent"],
        color="black",
    )

    ax.set_xlim(-0.5, 1.5)
    ax.set_xticks([0, 1], labels=CONDITIONS)

    ax.set_ylim(0, 1)
    ax.set_ylabel("Average CpG methylation", fontweight="bold")
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

    ax.spines[["top", "right"]].set_visible(False)

    ax.set_title(region, fontweight="bold", pad=12)

    fig.tight_layout()
    fig.savefig(
        f"{OUTPUT}/average_methylation-{region}.svg",
        bbox_inches="tight",
    )
    plt.close(fig)
