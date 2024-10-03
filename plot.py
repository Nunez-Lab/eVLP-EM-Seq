# %% Import

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import polars as pl

# %% Load data

raw_data = pl.read_csv("data/dss_output_OT.csv")

# %% Process data

data = (
    raw_data.with_columns(
        wt_mu=pl.col("mu1"),
        cd55off_mu=pl.col("mu2"),
        score=pl.when(pl.col("diff") > 0)  # when mu1 - mu2 > 0 (i.e. mu1 > mu2)
        # mu1 (WT) is greater, score should be negative > 0)
        .then(pl.col("pval").log10())
        # mu2 (CD55off) is greater, score should be positive
        .otherwise(pl.col("pval").log10().neg()),
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
        "wt_mu",
        "cd55off_mu",
        "score",
    )
    .sort(by=["chr_order", "pos"])
    .with_row_index()
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
        assert df[feature].min() > min(yticks)
        assert df[feature].max() < max(yticks)

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

        ax.axvline(
            lo,
            color="lightgray",
            linewidth=0.5,
        )

    left = df[index].min()
    right = df[index].max()

    ax.axvline(right, color="lightgray", linewidth=0.5)
    ax.set_xlim(left, right)

    if xlabel:
        ax.set_xlabel(xlabel)

    ax.set_xticks(
        xticks,
        labels=xticklabels,
        rotation=xtick_rotation,
    )

    ax.set_ylabel(
        ylabel if ylabel is not None else feature,
    )

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

    fig.savefig(f"graphs/controls_{condition}.png", dpi=300)
    plt.close(fig)

# %% Plot differential methylation

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

yticks = np.arange(-160, 161, 40)

for filter_expr, use_xticks, base_filename, sample in [
    (pl.col("chr_order") < 25, True, "differential_methylation_all", False),
    (CD55_SURROUNDING_EXPR, False, "differential_methylation_cd55", False),
]:
    df = data.sample(10_000) if sample else data
    fig, ax = manhattan(
        df.filter(filter_expr),
        by="chr",
        feature="score",
        two_sided=("CD55off", "WT"),
        highlight=None,
        yticks=yticks,
        yticklabels=abs(yticks),
        ylabel=r"$\bf{Significance}$ $\bf{score}$"
        + "\n"
        + r"$-\log_{10}(p$-value$)$",
        use_xticks=use_xticks,
    )

    ax.annotate(
        "CD55" + (" [SUBSAMPLED DATA]" if sample else ""),
        xy=(CD55_TSS_INDEX, 1.01),
        xytext=(CD55_TSS_INDEX, 1.12),
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

    suffix = "-SUBSAMPLED" if sample else ""

    fig.tight_layout()
    fig.savefig(f"graphs/{base_filename}{suffix}.png", dpi=300)
    plt.close(fig)
