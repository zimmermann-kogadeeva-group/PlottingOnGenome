from copy import deepcopy

from dna_features_viewer import BiopythonTranslator, GraphicRecord
from matplotlib.pyplot import Line2D


def shift_feature(feature, shift=0):
    """Helper function to shift a Biopython feature without changing the original"""
    new_feature = deepcopy(feature)
    new_feature.location = feature.location + shift
    return new_feature


def get_gene_coverage(gene_start, gene_end, region_start, region_end):
    if gene_start >= region_start and gene_end <= region_end:
        return 1.0
    elif gene_end < region_start or gene_start > region_end:
        return 0.0
    elif gene_start < region_start and gene_end > region_end:  # TODO: check this
        return (region_end - region_start + 1) / (gene_end - gene_start + 1)
    elif gene_start < region_start and region_start <= gene_end <= region_end:
        return (gene_end - region_start + 1) / (gene_end - gene_start + 1)
    elif gene_end > region_end and region_start <= gene_start <= region_end:
        return (region_end - gene_start + 1) / (gene_end - gene_start + 1)
    else:
        raise RuntimeError(
            f"Unknown gene case: "
            f"{gene_start=}, {gene_end=}, "
            f"insert_start={region_start}, insert_end={region_end}"
        )


def get_genes_graphic_record(
    genes,
    region_start,
    region_end,
    buffer=4000,
    col="#ebf3ed",
    feature_types=None,
    cmap=None,
):
    # Check args are of the correct type
    if isinstance(buffer, int):
        buffer = (buffer, buffer)
    assert isinstance(buffer, (tuple, list)) and len(buffer) == 2

    if feature_types is None:
        feature_types = {x.type for x in genes}

    if cmap is not None:
        for i, gene in enumerate(genes):
            genes[i].qualifiers["color"] = cmap(
                get_gene_coverage(
                    gene.location.start, gene.location.end, region_start, region_end
                )
            )
    conv = BiopythonTranslator()
    conv.default_feature_color = col
    features = [conv.translate_feature(x) for x in genes if x.type in feature_types]

    # Plot the genes and CDSes in the region of the mapped sequence
    record_hits = GraphicRecord(
        first_index=region_start - buffer[0],
        sequence_length=region_end - region_start + sum(buffer),
        features=features,
    )

    return record_hits


def fig_axvline(axes, value, ls="--", color="gray", **kwargs):
    fig = axes[0].get_figure()
    transFigure = fig.transFigure.inverted()

    coord1 = transFigure.transform(axes[0].transData.transform([value, 0]))
    coord2 = transFigure.transform(axes[1].transData.transform([value, 0]))

    line = Line2D(
        (coord1[0], coord2[0]),
        (coord1[1], coord2[1]),
        ls=ls,
        color=color,
        transform=fig.transFigure,
        zorder=-100,
    )
    fig.lines.append(line)
