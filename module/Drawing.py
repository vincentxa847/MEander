import os
import argparse
import logging
import polars as pl
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature
from dna_features_viewer import GraphicRecord

'''
Drawing module

INPUT  : gene / variant_type / variant_chr / variant_pos / variant_length / variant_strand / output_dir
OUTPUT : {variant_chr}_{variant_pos}_{variant_type}.png
'''

LOGGER = logging.getLogger(__name__)

# MEI is represented by a fixed-width placeholder on plot (also the default length if not provided)
MEI_PLACEHOLDER = 1

def strand(strand):
    """Convert '+'/'-' strand notation to DNA Features Viewer values."""
    if strand == "+":
        return 1
    elif strand == "-":
        return -1
    else:
        return 0

def extract_gff3_attribute(column, key):
    """Extract a GFF3 attribute such as gene_name or transcript_id."""
    return (
        pl.col(column)
        .str.split(";")
        .list.eval(
            pl.element()
            .filter(
                pl.element()
                .str.starts_with(f"{key}=")
            )
        )
        .list.first()
        .str.replace(
            f"{key}=",
            ""
        )
    )

def load_variants(
    path,
    melt_meinfo,
    melen_field,
    mestrand_field,
    default_length=MEI_PLACEHOLDER,
    default_strand=0,
):
    """Load and normalize the variant table."""

    variants = pl.read_csv(path)

    if melt_meinfo:
        if "SVLEN" in variants.columns:
            variant_length = (
                pl.col("SVLEN")
                .cast(pl.Int64, strict=False)
                .fill_null(default_length)
            )
        else:
            variant_length = pl.lit(default_length)

        if "MEINFO" in variants.columns:
            variant_strand = (
                pl.col("MEINFO")
                .cast(pl.String)
                .str.split(",")
                .list.get(3, null_on_oob=True)
                .fill_null(str(default_strand))
            )
        else:
            variant_strand = pl.lit(default_strand)

    else:
        if melen_field and melen_field in variants.columns:
            variant_length = (
                pl.col(melen_field)
                .cast(pl.Int64, strict=False)
                .fill_null(default_length)
            )
        else:
            variant_length = pl.lit(default_length)

        if mestrand_field and mestrand_field in variants.columns:
            variant_strand = (
                pl.col(mestrand_field)
                .fill_null(default_strand)
                .cast(pl.String)
            )
        else:
            variant_strand = pl.lit(default_strand)

    variants = variants.with_columns(
        variant_length.alias("variant_length"),
        variant_strand.alias("variant_strand"),
    )

    # Set gene from mapped pathogenic MEI or exclude additional information.
    variants = variants.with_columns(
        pl.when(pl.col("mapped_pathogenic_MEI") != ".")
        # from mapped pathogenic MEI
        .then(
            pl.col("mapped_pathogenic_MEI")
                .str.extract(r'Gene:"([^"]+)"', 1)
        )
        # exclude additional information
        .otherwise(
            pl.col("gene")
                .str.extract(r'^([^(]+)', 1)
        )
        .alias("gene")
    )

    return variants

def load_gff(path, genes):
    """Load MANE_Select annotations for genes in the variant table."""

    gff = (
        pl.scan_csv(
            path,
            separator="\t",
            comment_prefix="#",
            has_header=False,
            new_columns=[
                "chrom",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attribute",
            ],
        )
        .filter( # subset to MANE_Select to reduce load
            pl.col("attribute").str.contains("MANE_Select")
        )
        .collect()
    )

    gff = gff.with_columns(
        extract_gff3_attribute(
            "attribute",
            "gene_name",
        ).alias("gene_name"),

        extract_gff3_attribute(
            "attribute",
            "transcript_id",
        ).alias("transcript_id"),

        extract_gff3_attribute(
            "attribute",
            "exon_number",
        ).alias("exon_number"),

        extract_gff3_attribute(
            "attribute",
            "tag",
        ).alias("tag"),
    )

    return gff.filter(
        pl.col("gene_name").is_in(genes)
    )

def make_features(
    gene_start,
    exons,
    label_exon=None,
    introns=None,
    label_intron=None,
):
    """Convert exon/intron coordinates into DNA Features Viewer features."""

    features = []

    # Exons
    for row in exons.iter_rows(named=True):

        label = (
            f"Exon {row['exon_number']}"
            if row["exon_number"] == label_exon
            else None
        )

        features.append(
            GraphicFeature(
                start=row["start"] - gene_start,
                end=row["end"] - gene_start,
                strand=strand(row["strand"]),
                label=label,
                color="skyblue",
            )
        )

    # Intron
    if label_intron is not None and introns.height > 0:

        selected = introns.filter(
            pl.col("intron_number") == label_intron
        )

        if selected.height > 0:
            row = selected.row(0, named=True)

            # Small feature as the label anchor.
            intron_mid = (
                row["start"] + row["end"]
            ) // 2

            features.append(
                GraphicFeature(
                    start=intron_mid - gene_start,
                    end=intron_mid - gene_start,
                    label=f"Intron {row['intron_number']}",
                    thickness=0,
                )
            )

    return features

def build_introns(exons, tx_strand):
    """Build introns between transcript-strand-ordered exons."""

    introns = []

    for i in range(exons.height - 1):

        if tx_strand == "+":
            intron_start = exons["end"][i]
            intron_end = exons["start"][i + 1]
        else:
            intron_start = exons["end"][i + 1]
            intron_end = exons["start"][i]

        introns.append(
            {
                "start": min(intron_start, intron_end),
                "end": max(intron_start, intron_end),
                "intron_number": i + 1,
            }
        )

    if not introns:
        return pl.DataFrame(
            schema={
                "start": pl.Int64,
                "end": pl.Int64,
                "intron_number": pl.Int64,
            }
        )

    return pl.DataFrame(introns)


def get_variant_context(exons, introns, variant_pos, tx_strand):
    """
    Determine whether the variant is exonic or intronic and
    identify the appropriate exon/intron label.
    """

    variant_exon = exons.filter(
        (pl.col("start") <= variant_pos)
        & (pl.col("end") >= variant_pos)
    )

    if variant_exon.height > 0:
        exon_number = variant_exon["exon_number"][0]

        return {
            "in_exon": True,
            "variant_exon": variant_exon,
            "exon_label": exon_number,
            "intron_label": None,
            "overview_label": f"Exon {exon_number}",
        }

    affected_introns = introns.filter(
        (pl.col("start") <= variant_pos)
        & (pl.col("end") >= variant_pos)
    )

    if affected_introns.height == 0:
        return None

    affected_intron = affected_introns.row(0, named=True)
    intron_number = affected_intron["intron_number"]

    return {
        "in_exon": False,
        "variant_exon": variant_exon,
        "exon_label": None,
        "intron_label": intron_number,
        "overview_label": f"Intron {intron_number}",
    }


def make_zoom_record(
    exons,
    introns,
    context,
    variant_pos,
    variant_type,
    variant_length,
    variant_strand,
    gene_start,
    zoom_padding=200,
):
    """Create the zoomed GraphicRecord."""

    zoom_features = []

    if context["in_exon"]:

        # NOTE: MEI overlap with exon feature and MEI will be push up
        affected_exon = context["variant_exon"].row(0, named=True)

        exon_start = affected_exon["start"]
        exon_end = affected_exon["end"]

        zoom_start = exon_start - zoom_padding
        zoom_end = exon_end + zoom_padding

        zoom_features.append(
            GraphicFeature(
                start=exon_start - zoom_start,
                end=exon_end - zoom_start,
                strand=strand(affected_exon["strand"]),
                label=f"Exon {affected_exon['exon_number']}",
                color="skyblue",
            )
        )

    else:

        affected_intron = introns.filter(
            pl.col("intron_number")
            == context["intron_label"]
        ).row(0, named=True)

        intron_start = affected_intron["start"]
        intron_end = affected_intron["end"]

        zoom_start = max(
            intron_start - zoom_padding,
            gene_start,
        )

        zoom_end = intron_end + zoom_padding

        neighboring_exons = exons.filter(
            (pl.col("end") >= zoom_start)
            & (pl.col("start") <= zoom_end)
        )

        for row in neighboring_exons.iter_rows(named=True):

            zoom_features.append(
                GraphicFeature(
                    start=max(
                        row["start"],
                        zoom_start,
                    ) - zoom_start,
                    end=min(
                        row["end"],
                        zoom_end,
                    ) - zoom_start,
                    strand=strand(row["strand"]),
                    label=f"Exon {row['exon_number']}",
                    color="skyblue",
                )
            )

    mei_label = (
        f"{variant_type}"
        if variant_length == MEI_PLACEHOLDER
        else f"{variant_type} ({variant_length} bp)"
    )

    zoom_features.append(
        GraphicFeature(
            start=variant_pos - zoom_start,
            end=variant_pos - zoom_start + MEI_PLACEHOLDER,
            strand=strand(variant_strand),
            label=mei_label,
            color="red",
        )
    )

    return GraphicRecord(
        sequence_length=zoom_end - zoom_start,
        features=zoom_features,
    )


def plot_variant(
    variant,
    gene_gff,
    output_dir,
):
    """Generate one plot for one variant."""

    gene = variant["gene"]
    variant_type = variant["ALT"]
    variant_chr = variant["CHROM"]
    variant_pos = int(variant["POS"])
    variant_length = int(variant["variant_length"])
    variant_strand = variant["variant_strand"]

    exons = gene_gff.filter(
        pl.col("feature") == "exon"
    )

    if exons.height == 0:
        LOGGER.warning(
            "No exon annotation found for %s; skipping.",
            gene,
        )
        return

    tx_strand = exons["strand"][0]

    # Transcript-strand-aware exon order
    if tx_strand == "+":
        exons = exons.sort("start")
    else:
        exons = exons.sort(
            "start",
            descending=True,
        )

    introns = build_introns(
        exons,
        tx_strand,
    )

    gene_start = exons["start"].min()
    gene_end = exons["end"].max()

    context = get_variant_context(
        exons,
        introns,
        variant_pos,
        tx_strand,
    )

    if context is None:
        LOGGER.warning(
            "Variant %s:%s does not overlap an exon/intron of %s; skipping.",
            variant_chr,
            variant_pos,
            gene,
        )
        return

    # -----------------------------
    # Overview
    # -----------------------------

    gene_features = make_features(
        gene_start=gene_start,
        exons=exons,
        label_exon=context["exon_label"],
        introns=introns,
        label_intron=context["intron_label"],
    )

    gene_record = GraphicRecord(
        sequence_length=gene_end - gene_start,
        features=gene_features,
    )

    # -----------------------------
    # Zoom
    # -----------------------------

    zoom_record = make_zoom_record(
        exons=exons,
        introns=introns,
        context=context,
        variant_pos=variant_pos,
        variant_type=variant_type,
        variant_length=variant_length,
        variant_strand=variant_strand,
        gene_start=gene_start,
    )

    # -----------------------------
    # Plot
    # -----------------------------

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(15, 5),
        gridspec_kw={
            "height_ratios": [1, 2]
        },
    )

    fig.subplots_adjust(
        hspace=0.8
    )

    # Overview
    gene_record.plot(
        ax=axes[0],
        figure_width=15,
    )

    axes[0].set_title(
        f"{gene}"
    )

    # Zoom
    zoom_record.plot(
        ax=axes[1],
        figure_width=15,
    )

    axes[1].set_title(
        context["overview_label"]
    )

    # -----------------------------
    # Export
    # -----------------------------

    filename = (
        f"{variant_chr}_"
        f"{variant_pos}_"
        f"{variant_type}.png"
    )

    output_path = os.path.join(
        output_dir,
        filename,
    )

    fig.savefig(
        output_path,
        dpi=150,
        bbox_inches="tight",
    )

    plt.close(fig)

    LOGGER.debug("Plot: %s", output_path)


def main(argv=None):

    parser = argparse.ArgumentParser(
            description="Draw gene/transcript structures and MEI variant locations."
        )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Annotated variant CSV file to plot.",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Directory for output PNG files.",
    )

    parser.add_argument(
        "-g",
        "--gff",
        required=True,
        help="Gencode GFF3 annotation file.",
    )

    parser.add_argument(
        "--melt-meinfo",
        action="store_true",
        help="Handle MEI length and strand using MELT MEINFO fields.",
    )

    parser.add_argument(
        "-melen",
        "--melen-field",
        default=None,
        help="Column of the MEI length.",
    )
    
    parser.add_argument(
        "-mestrand",
        "--mestrand-field",
        default=None,
        help="Column of the MEI strand.",
    )

    args = parser.parse_args(argv)


    os.makedirs(
        args.output_dir,
        exist_ok=True,
    )

    ## Load variants
    variants = load_variants(
        path=args.input,
        melt_meinfo=args.melt_meinfo,
        melen_field=args.melen_field,
        mestrand_field=args.mestrand_field,
    )

    ## Load GFF
    genes = variants["gene"].drop_nulls().unique()

    gff = load_gff(
        path=args.gff,
        genes=genes,
    )

    ## Plot
    for variant in variants.iter_rows(named=True):

        gene = variant["gene"]

        if gene is None:
            LOGGER.warning(
                "variant has no gene; skipping."
            )
            continue

        gene_gff = gff.filter(
            pl.col("gene_name") == gene
        )

        if gene_gff.height == 0:
            LOGGER.warning("Gene %s not found in GFF; skipping.", gene)
            continue

        plot_variant(
            variant=variant,
            gene_gff=gene_gff,
            output_dir=args.output_dir,
        )


if __name__ == "__main__":
    main()
