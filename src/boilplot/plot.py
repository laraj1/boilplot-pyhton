import numpy as np
import pandas as pd

from plotnine import (
    aes,
    facet_wrap,
    geom_hline,
    geom_jitter,
    geom_rect,
    geom_text,
    geom_violin,
    ggplot,
    labs,
    scale_fill_manual,
    scale_x_continuous,
    scale_y_continuous,
    theme_classic,
)

from mizani.breaks import breaks_extended
from mizani.palettes import hue_pal


# General utilities

def nonnegative_breaks(limits):
    """Generate y-axis breaks restricted to non-negative values."""
    breaks = breaks_extended(n=5)(limits)
    return [b for b in breaks if b >= 0]


# Input validation

def _validate_inputs(
    adata,
    genes,
    category_col,
    category_order,
    palette,
    min_nonzero,
    layer,
    use_raw,
):
    """Validate arguments supplied to boil_plot."""

    if layer is not None and use_raw:
        raise ValueError("layer and use_raw cannot both be specified.")

    if layer is not None and layer not in adata.layers:
        raise ValueError(
            f"Layer '{layer}' was not found in adata.layers."
        )

    if use_raw and adata.raw is None:
        raise ValueError(
            "use_raw=True was requested, but adata.raw is None."
        )

    if not isinstance(min_nonzero, (int, np.integer)) or min_nonzero < 1:
        raise ValueError("min_nonzero must be a positive integer.")

    if category_col not in adata.obs.columns:
        raise ValueError(f"'{category_col}' was not found in adata.obs.")

    if isinstance(genes, str):
        genes = [genes]
    else:
        genes = list(genes)

    if len(genes) == 0:
        raise ValueError("genes must contain at least one gene.")

    genes = list(dict.fromkeys(genes))

    if use_raw:
        available_genes = adata.raw.var_names
        source = "adata.raw.var_names"
    else:
        available_genes = adata.var_names
        source = "adata.var_names"

    missing_genes = [
        gene
        for gene in genes
        if gene not in available_genes
    ]

    if missing_genes:
        raise ValueError(
            f"The following genes were not found in "
            f"{source}: {missing_genes}"
        )

    if category_order is not None:
        category_order = list(category_order)

        if len(category_order) == 0:
            raise ValueError("category_order cannot be empty.")

        available_categories = set(
            adata.obs[category_col].dropna().unique()
        )

        missing_categories = [
            category
            for category in category_order
            if category not in available_categories
        ]

        if missing_categories:
            raise ValueError(
                "The following categories in category_order "
                f"were not found in adata.obs['{category_col}']: "
                f"{missing_categories}"
            )

    if palette is not None and not isinstance(palette, dict):
        raise ValueError(
            "palette must be either a dictionary or None."
        )

    return genes


# Category handling

def _get_categories(
    adata,
    category_col,
    category_order=None,
):
    """Determine the categories and their plotting order."""

    values = adata.obs[category_col]

    if category_order is not None:
        return list(category_order)

    if (
        isinstance(values.dtype, pd.CategoricalDtype)
        and values.cat.ordered
    ):
        observed = set(values.dropna().unique())

        return [
            category
            for category in values.cat.categories
            if category in observed
        ]

    return list(values.dropna().unique())


# Palette handling

def _get_palette(
    palette,
    categories,
):
    """Convert palette specification into a category -> colour dict."""

    if isinstance(palette, dict):
        missing = [
            category
            for category in categories
            if category not in palette
        ]

        if missing:
            raise ValueError(
                "The supplied palette is missing colours for "
                f"the following categories: {missing}"
            )

        return {
            category: palette[category]
            for category in categories
        }

    colours = hue_pal()(len(categories))

    return dict(zip(categories, colours))


# Expression extraction

def _get_expression(
    adata,
    genes,
    layer=None,
    use_raw=False,
):
    """Extract expression for the requested genes."""

    if use_raw:
        expr = adata.raw[:, genes].X
    elif layer is not None:
        expr = adata[:, genes].layers[layer]
    else:
        expr = adata[:, genes].X

    if hasattr(expr, "toarray"):
        expr = expr.toarray()

    return np.asarray(expr)


# Plot data preparation

def _prepare_plot_data(
    adata,
    genes,
    category_col,
    categories,
    expression,
    min_nonzero,
):
    """Prepare long-format data and nonzero-proportion information."""

    df = pd.DataFrame(
        expression,
        columns=genes,
        index=adata.obs_names,
    )

    df[category_col] = adata.obs[category_col].values

    df = df[
        df[category_col].isin(categories)
    ].copy()

    if isinstance(df[category_col].dtype, pd.CategoricalDtype):
        df[category_col] = pd.Categorical(
            df[category_col],
            categories=categories,
            ordered=True,
        )

    expression_long = df.melt(
        id_vars=category_col,
        value_vars=genes,
        var_name="gene",
        value_name="expression",
    )

    # Preserve the order in which genes were supplied by the user
    expression_long["gene"] = pd.Categorical(
        expression_long["gene"],
        categories=genes,
        ordered=True,
    )

    category_positions = {
        category: float(i)
        for i, category in enumerate(categories)
    }

    expression_long["x"] = (
        expression_long[category_col]
        .map(category_positions)
        .astype(float)
    )

    nonzero_counts = (
        expression_long
        .groupby(
            ["gene", category_col],
            observed=True,
        )["expression"]
        .apply(
            lambda x: (x > 0).sum()
        )
        .reset_index(
            name="nonzero_count"
        )
    )

    expression_long = expression_long.merge(
        nonzero_counts,
        on=["gene", category_col],
        how="left",
    )

    detection_data = (
        expression_long
        .groupby(
            ["gene", category_col],
            observed=True,
        )
        .agg(
            nonzero_prop=(
                "expression",
                lambda x: (x > 0).mean()
            )
        )
        .reset_index()
    )

    detection_data["x"] = (
        detection_data[category_col]
        .map(category_positions)
        .astype(float)
    )

    max_expression = (
        expression_long
        .groupby(
            "gene",
            observed=True,
        )["expression"]
        .max()
        .reset_index(
            name="max_expression"
        )
    )

    detection_data = detection_data.merge(
        max_expression,
        on="gene",
        how="left",
    )

    bar_width = 0.65

    detection_data["bar_height"] = (
        detection_data["max_expression"] * 0.15
    )

    detection_data["bar_ymin"] = (
        -detection_data["bar_height"]
    )

    detection_data["bar_ymax"] = (
        -detection_data["max_expression"] * 0.05
    )

    detection_data["bar_xmin"] = (
        detection_data["x"] - bar_width / 2
    )

    detection_data["bar_xmax"] = (
        detection_data["x"] + bar_width / 2
    )

    detection_data["detection_xmax"] = (
        detection_data["bar_xmin"]
        + bar_width * detection_data["nonzero_prop"]
    )

    detection_data["label"] = (
        (
            detection_data["nonzero_prop"] * 100
        )
        .round()
        .astype(int)
        .astype(str)
        + "%"
    )

    detection_data["label_y"] = (
        -detection_data["bar_height"] * 1.5
    )

    strip_data = expression_long[
        expression_long["expression"] > 0
    ]

    violin_data = expression_long[
        (expression_long["expression"] > 0)
        &
        (expression_long["nonzero_count"] >= min_nonzero)
    ]

    return (
        expression_long,
        strip_data,
        violin_data,
        detection_data,
        category_positions,
    )


# Plot construction

def _build_plot(
    expression_long,
    strip_data,
    violin_data,
    detection_data,
    category_positions,
    category_col,
    palette,
    title,
    point_size,
    point_alpha,
    jitter_width,
):
    """Construct the plotnine plot."""

    p = (
        ggplot(
            expression_long,
            aes(
                x="x",
                y="expression",
                fill=category_col,
            )
        )
        + geom_violin(
            data=violin_data,
            trim=True,
            scale="width",
            width=0.9,
        )
        + geom_jitter(
            data=strip_data,
            alpha=point_alpha,
            size=point_size,
            width=jitter_width,
        )
        + geom_rect(
            data=detection_data,
            mapping=aes(
                xmin="bar_xmin",
                xmax="bar_xmax",
                ymin="bar_ymin",
                ymax="bar_ymax",
            ),
            fill="#bdbdbd",
            inherit_aes=False,
        )
        + geom_rect(
            data=detection_data,
            mapping=aes(
                xmin="bar_xmin",
                xmax="detection_xmax",
                ymin="bar_ymin",
                ymax="bar_ymax",
                fill=category_col,
            ),
            inherit_aes=False,
        )
        + geom_text(
            data=detection_data,
            mapping=aes(
                x="x",
                y="label_y",
                label="label",
            ),
            size=7,
            color="#333333",
            inherit_aes=False,
        )
        + geom_hline(
            yintercept=0,
            color="#444444",
            size=0.3,
        )
        + facet_wrap(
            "~gene",
            scales="free_y",
            ncol=1,
        )
        + labs(
            title=title,
            fill=category_col,
            x=category_col,
            y="Expression",
        )
        + scale_fill_manual(
            values=palette
        )
        + scale_x_continuous(
            breaks=list(category_positions.values()),
            labels=list(category_positions.keys()),
            expand=(0.1, 0.1),
        )
        + scale_y_continuous(
            breaks=nonnegative_breaks,
        )
        + theme_classic()
    )

    return p


# Public API

def boil_plot(
    adata,
    genes,
    category_col,
    palette=None,
    category_order=None,
    min_nonzero=10,
    layer=None,
    use_raw=False,
    title=None,
    point_size=0.7,
    point_alpha=0.6,
    jitter_width=0.1,
):
    """
    Plot single-cell expression distributions for one or more genes.

    Only non-zero expression values are shown as a violin + jitter plot.
    The violin is only shown for gene/category combinations with at
    least `min_nonzero` non-zero observations.

    The proportion of non-zero observations is shown as a bar below
    zero.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing expression data.

    genes : str or list of str
        Gene or genes to plot.

    category_col : str
        Column in adata.obs containing the categories.

    palette : dict or None, default=None
        Colour specification.

        A dictionary should map category names to colours.

        If None, a categorical palette is generated automatically.

    category_order : list or None, default=None
        Optional explicit order of categories on the x-axis.

    min_nonzero : int, default=10
        Minimum number of non-zero observations required for a
        violin plot to be shown for a gene/category combination.

    layer : str or None, default=None
        AnnData layer from which to extract expression.

        If None, adata.X is used.

    use_raw : bool, default=False
        If True, expression is extracted from adata.raw.

    title : str or None, default=None
        Optional plot title.

    point_size : float, default=0.7
        Size of strip plot points.

    point_alpha : float, default=0.6
        Transparency of strip plot points.

    jitter_width : float, default=0.1
        Horizontal jitter applied to strip plot points.

    Returns
    -------
    plotnine.ggplot
        The resulting plot.
    """

    genes = _validate_inputs(
        adata=adata,
        genes=genes,
        category_col=category_col,
        category_order=category_order,
        palette=palette,
        min_nonzero=min_nonzero,
        layer=layer,
        use_raw=use_raw,
    )

    categories = _get_categories(
        adata=adata,
        category_col=category_col,
        category_order=category_order,
    )

    resolved_palette = _get_palette(
        palette=palette,
        categories=categories,
    )

    expression = _get_expression(
        adata=adata,
        genes=genes,
        layer=layer,
        use_raw=use_raw,
    )

    (
        expression_long,
        strip_data,
        violin_data,
        detection_data,
        category_positions,
    ) = _prepare_plot_data(
        adata=adata,
        genes=genes,
        category_col=category_col,
        categories=categories,
        expression=expression,
        min_nonzero=min_nonzero,
    )

    return _build_plot(
        expression_long=expression_long,
        strip_data=strip_data,
        violin_data=violin_data,
        detection_data=detection_data,
        category_positions=category_positions,
        category_col=category_col,
        palette=resolved_palette,
        title=title,
        point_size=point_size,
        point_alpha=point_alpha,
        jitter_width=jitter_width,
    )
