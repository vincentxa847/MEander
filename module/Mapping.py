import argparse
import os
import re
import tempfile
import logging
from collections import defaultdict
from contextlib import contextmanager
import polars as pl
import pybedtools
from cyvcf2 import VCF

from other_database import prepare_vcf
from module.QC import passes_qc # ModuleNotFoundError: No module named 'module.QC'

'''
Apply QC and annotation modules

Mapping module
- Genomic features of MANE Select transcripts (CDS, UTR, junction)
- Disease-focused MEIs within the defined overlap tolerance
- ClinGen information and gene pLI for MEI-overlapping genes
- gnomAD allele frequency

INPUT  : VCF or a list with VCFs (must be GRCh38)
OUTPUT : annotated table
'''

LOGGER = logging.getLogger(__name__)

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Map and annotate MEI VCF files."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input VCF or .list.txt file"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output CSV"
    )

    parser.add_argument(
        "--position-to-map",
        choices=[
                "start",
                "end",
                ],
        default="start",
        help="Position in disease-associated MEI catalog (TSD start or end) to use."
    )

    parser.add_argument(
        "--mane",
        required= True,
        help="BED file of MANE Select Transcript features."
    )

    parser.add_argument(
        "--pathogenic",
        required= True,
        help="Disease focused MEI catalog."
    )

    parser.add_argument(
        "--clingen",
        required= True,
        help="ClinGen Gene-Disease Validity tables."
    )

    parser.add_argument(
        "--gnomad-constraint",
        required= True,
        help="gnomad constraint metrics."
    )

    parser.add_argument(
        "--control-dir",
        required= True,
        help="Directory of gnomad MEIs."
    )

    parser.add_argument(
        "--mode",
        nargs="+",
        choices=[
            "CDS",
            "junction_donor",
            "junction_acceptor",
            "three_prime_UTR",
            "five_prime_UTR",
        ],
        default=["CDS"],
        help="Keep MEI in selected genomic features."
    )

    parser.add_argument(
        "--clingen-only",
        action="store_true",
        help="Keep MEI in genes with ClinGen Gene-Disease Clinical Validity",
    )

    parser.add_argument(
        "--qc",
        choices=[
            "MELT",
            "MEGAnE",
            "xTea",
        ],
        default=None,
        help="QC to run.",
    )

    parser.add_argument(
        "--siteonly",
        action="store_true",
        help="When enabled, sample-level information is not collected.",
    )

    parser.add_argument(
        "--window",
        type=int,
        default=20,
        help="Overlap tolerance"
    )

    parser.add_argument(
        "--alt",
        default="ALT",
        help="VCF field of the MEI."
    )

    parser.add_argument(
        "--clingen-classification",
        nargs="+",
        choices=[
            "Definitive",
            "Strong",
            "Moderate",
            "Limited",
            "No Known Disease Relationship",
            "Disputed",
            "Refuted",
        ],
        default=["Definitive", "Strong"],
        help="ClinGen classification."
    )

    return parser.parse_args(argv)


# control cohort
# NOTE: some polymorphic MEIs not being classified as any ME family (i.e. '<INS>' not '<INS:ME:ALU>'), therefore none in gnomAD
control_vcfs = ["gnomad.v4.1.sv.sites.INS.ALU.PASS.vcf.gz", "gnomad.v4.1.sv.sites.INS.LINE1.PASS.vcf.gz",
                "gnomad.v4.1.sv.sites.INS.SVA.PASS.vcf.gz", "gnomad.v4.1.sv.sites.INS.HERVK.PASS.vcf.gz"]
control_af   = "AF"
POP_AF_FIELDS = {
    "AF_afr": "afr",
    "AF_ami": "ami",
    "AF_amr": "amr",
    "AF_asj": "asj",
    "AF_eas": "eas",
    "AF_fin": "fin",
    "AF_mid": "mid",
    "AF_nfe": "nfe",
    "AF_rmi": "rmi",
    "AF_sas": "sas",
}

# harmonize MEI representation
TE_ALIASES = {
    "Alu"   : ["INS:ME:ALU", "alu_insertion", "Alu"],
    "LINE-1": ["INS:ME:LINE1", "line1 insertion", "L1"],
    "SVA"   : ["INS:ME:SVA", "sva insertion", "SVA"],
}

@contextmanager
def redirect_htslib_stderr():
    """
    Capture cyvcf2/HTSlib messages written directly to OS-level stderr (fd2)
    """
    with tempfile.TemporaryFile(mode="w+") as f:

        old_stderr = os.dup(2)

        try:
            os.dup2(f.fileno(), 2)
            yield

        # regardless of whether VCF() succeeds or raises an exception, restores stderr
        finally:
            os.dup2(old_stderr, 2)
            os.close(old_stderr)

        f.seek(0)
        msg = f.read()

    if msg.strip():
        for line in msg.splitlines():
            LOGGER.warning("cyvcf2/HTSlib: %s", line)


def open_vcf(path):
    with redirect_htslib_stderr():
        vcf = VCF(path)

    return vcf

def harmonize_te(value):
    value = str(value)
    for te, aliases in TE_ALIASES.items():
        if any(alias.lower() in str(value).lower() for alias in aliases):
            return te
    return value

# feature rank
feature_order = {
    "CDS"               : 1,
    "junction_donor"    : 2,
    "junction_acceptor" : 2,
    "three_prime_UTR"   : 3,
    "five_prime_UTR"    : 4,
}

# if multiple features exist (e.g. CDS;junction), use the highest priority
def feature_rank(x):
    if x is None or x == "":
        return 99

    ranks = [
        feature_order.get(f, 99)
        for f in x.split(";")
    ]

    return min(ranks)


def mapping(
    chrom,
    pos,
    alt,
    window,
    database,
    mode
):
    hits = []

    for x in database.get(chrom, []):

        if harmonize_te(x["alt"]) != alt:
            continue
        
        if abs(x["pos"] - pos) <= window:

          if mode == "gnomad":
                hits.append(
                    f"{chrom}:{x['pos']}|{x['alt']}|{x['af']}|"
                    f"MAF:{x['grpmax_af']}|MAF_POP:{x['grpmax_pop']}"
                )

          elif mode == "pathogenic":
                hits.append(
                    x["annotation"]
                )

    return ";".join(hits) if hits else "."


def mane_intersection(mane, vcf, mode, alt_field):
    ''' 
    Map with MANE genomic feature 
        1. convert VCF to BED
        2. bedtools intersection
    '''

    # 1. convert VCF to BED
    vcf_records = []

    for v in vcf:

        if alt_field.startswith("INFO/"):
            info_key = alt_field.split("/", 1)[1]
            alt_values = [v.INFO.get(info_key)]
        else:
            alt_values = getattr(v, alt_field)

        for alt in alt_values:
            if alt is None:
                continue

            vcf_records.append([
                v.CHROM,
                v.POS - 1,
                v.POS,
                v.ID or ".",
                v.REF,
                str(alt)
            ])

    vcf_bed = pybedtools.BedTool(vcf_records)

    # 2. bedtools intersection
    mane_bed = pybedtools.BedTool(mane)

    mane_hits = vcf_bed.intersect(
        mane_bed,
        wa=True,
        wb=True
    )

    in_mane_annotation = defaultdict(list)

    for x in mane_hits:

        f = x.fields

        key = (
            str(f[0]),            # chr
            str(int(f[1]) + 1),   # POS (1-based)
            str(f[4]),            # REF
            str(f[5])             # ALT or SVTYPE
        )

        tid     = f[9]
        strand  = f[10]
        feature = f[11]
        gene    = f[12]
        exon    = f[13]

        # NOTE: pathogenic MEI in intron or not in mode should fetch annotation from 'mapped_pathogenic_MEI'  
        if feature not in mode:
            continue

        in_mane_annotation[key].append(
            {
                "gene"    : f"{gene}({tid},{strand})",
                "feature" : feature,
                "exon"    : exon
            }
        )

    return in_mane_annotation


def annotate_vcf(
    vcf,
    in_mane_annotation,
    alt_field,
    source,
    qc,
    clingen_only,
    siteonly,
    window,
    pathogenic_mei_index,
    control_af_index,
    clingen_dict,
    gnomad_constraint_dict,
):
    ''' 
    Main functionality
        1. apply QC module
        2. annotate MANE Select Transcript mapped
        3. annotate pathogenic MEI mapped
        4. annotate allele frequency in control cohort (gnomad)
        5. annotate clingen information
        6. annotate gnomAD constraint of transcript located

    Returns
        all_variants_of_this_file : list Variants that passed QC and downstream filtering.
        qc_failed_count : int Number of VCF records that failed QC.
    
    '''

    all_variants_of_this_file = []
    qc_failed_count = 0

    for v in vcf:

        # --- QC
        if qc is not None:
            if not passes_qc(v, qc):
                qc_failed_count += 1
                continue

        if alt_field.startswith("INFO/"):
            info_key = alt_field.split("/", 1)[1]
            alt_values = [v.INFO.get(info_key)]
        else:
            alt_values = getattr(v, alt_field)

        for alt in alt_values:

            key = (
                str(v.CHROM),
                str(v.POS),
                str(v.REF),
                str(alt)
            )

            alt = harmonize_te(alt)

            row = {
                "source": source,
                "CHROM" : v.CHROM,
                "POS"   : v.POS,
                "ID"    : v.ID,
                "REF"   : v.REF,
                "ALT"   : alt,
                "FILTER": v.FILTER
            }

            # --- INFO
            for info_key, value in v.INFO:

                if value is None:
                    row[info_key] = "."
                elif isinstance(value, (list, tuple)):
                    row[info_key] = ",".join(map(str, value))
                else:
                    row[info_key] = value

            # --- MANE

            if key in in_mane_annotation:

                info = in_mane_annotation[key]

                for col in ["gene", "feature", "exon"]:
                    row[col] = ";".join(
                        sorted({
                            x[col]
                            for x in info
                        })
                    )

            else:
                row["feature"] = "."
                row["gene"]    = "."
                row["exon"]    = "."


            # --- pathogenic MEI

            row["mapped_pathogenic_MEI"] = mapping(
                v.CHROM,
                v.POS,
                alt,
                window,
                pathogenic_mei_index,
                mode="pathogenic"
            )

            # fallback of gene annotation
            if row["gene"] == "." and row["mapped_pathogenic_MEI"] != ".":

                gene = re.findall(
                    r'Gene:"([^"]+)"',
                    row["mapped_pathogenic_MEI"]
                )

                transcript = re.findall(
                    r'Transcript:"(ENST[^"/\s]+)',
                    row["mapped_pathogenic_MEI"]
                )

                strand = re.findall(
                    r'Transcript:"[^"]+\(([-+])\)"',
                    row["mapped_pathogenic_MEI"]
                )

                if gene and transcript and strand:
                    row["gene"] = f"{gene[0]}({transcript[0]}, {strand[0]})"

            # --- gnomad AF

            row["gnomad_AF"] = mapping(
                v.CHROM,
                v.POS,
                alt,
                window,
                control_af_index,
                mode="gnomad"
            )

            # --- clingen annotation

            gene_symbol = (
                row["gene"]
                .split("(")[0]
            )

            row["clingen annotation"] = clingen_dict.get(
                gene_symbol,
                "."
            )

            # --- gnomAD constraint

            match = re.search(
                r"\((ENST[^,]+)",
                row["gene"]
            )

            if match:

                tid = match.group(1).split(".")[0]
                lof = gnomad_constraint_dict.get(tid)

                # lof.pLI is None when tid does not exist as key, missing value in gnomad constraint encodes as 'NA' 
                row["lof.pLI"] = None if lof in (None, "NA") else round(float(lof), 4)

            else:
                row["lof.pLI"] = None

            # --- Sample information

            if not siteonly:
                gt = v.genotypes

                samples = []
                for i, sample_name in enumerate(vcf.samples):
                    gt_value = f"{gt[i][0]}/{gt[i][1]}"

                    # Keep only GT and samples carrying the allele
                    if gt_value in {"0/0", "./.", "-1/-1"}:
                        continue

                    sample_value = gt_value

                    samples.append(f"{sample_name}:{sample_value}")

                row["sample"] = ";".join(samples) if samples else "."

            # --- MEI to keep
            # 1. mapped pathogenic MEI
            # 2. in MANE features in mode
            # 3. ClinGen filtering
            keep = (
                row["mapped_pathogenic_MEI"] != "."
                or (
                    key in in_mane_annotation
                    and (
                        not clingen_only
                        or row["clingen annotation"] != "."
                    )
                )
            )

            if keep:
                all_variants_of_this_file.append(row)

    return all_variants_of_this_file, qc_failed_count

# ---------------------------------------------------------
# Main
# single cohort VCF : one variant in two samples get one rows
# a list of VCFs    : one variant in two samples get two rows
# ---------------------------------------------------------

def main(argv=None):

    args = parse_args(argv)

    LOGGER.debug("Arguments: %s", vars(args))

    mode           = args.mode
    siteonly       = args.siteonly
    window         = args.window
    patho_pos_2map = args.position_to_map
    clingen_only   = args.clingen_only

    # ---------------------------------------------------------
    # Input
    # ---------------------------------------------------------

    if args.input.endswith(".list.txt"):
        with open(args.input) as f:
            vcf_files = [
                line.strip()
                for line in f
                if line.strip()
            ]
    else:
        vcf_files = [args.input]

    LOGGER.debug("Input VCFs: %d", len(vcf_files))

    if siteonly:
        LOGGER.info("Site-only mode enabled: sample-level information will not be collected")

    # ---------------------------------------------------------
    # pathogenic_mei
    # ---------------------------------------------------------
    LOGGER.debug("Loading pathogenic MEI catalog: %s", args.pathogenic)

    pathogenic_mei = pl.read_excel(
        args.pathogenic,
        sheet_name="Insertion"
    )

    # NOTE: change here if column name change
    pathogenic_group      = "Group"
    pathogenic_start      = "Position (GRCh38)"
    pathogenic_end        = "END (GRCh38)"
    pathogenic_gene       = "Gene"
    pathogenic_transcript = "MANE Select transcript (orientation)"
    pathogenic_feature    = "Genomic Feature"
    pathogenic_te         = "Insertion"
    pathogenic_subfamily  = "Subfamily"
    pathogenic_strand     = "Strand"
    pathogenic_size       = "Size (bp)"
    pathogenic_disease    = "Phenotype (MIM number or MedGen UID)"
    pathogenic_moi        = "Mode of Inheritance"
    pathogenic_pmid       = "PMID or ClinVar ID"

    if patho_pos_2map == "start":
        pathogenic_pos = pathogenic_start
    elif patho_pos_2map == "end":
        pathogenic_pos = pathogenic_end

    # record with valid position
    pattern = r"^chr[^:]+:\d+$"
    pathogenic_mei = pathogenic_mei.with_columns(
        pl.col(pathogenic_pos)
        .cast(pl.String)
        .str.strip_chars()
        .alias("_position")
    )

    valid_position = (
        pl.col("_position")
        .str.contains(pattern)
    )

    dropped = pathogenic_mei.filter(
        ~valid_position
    ).height

    if dropped > 0:
        LOGGER.debug("Dropped %s pathogenic MEI records with invalid Position format", dropped)

    pathogenic_mei = (
        pathogenic_mei
        .filter(valid_position)
        .with_columns(
            [
                pl.col("_position")
                .str.split_exact(":", 1)
                .struct.field("field_0")
                .alias("CHROM"),

                pl.col("_position")
                .str.split_exact(":", 1)
                .struct.field("field_1")
                .cast(pl.Int64)
                .alias("POS"),
            ]
        )
        .drop("_position")
    )

    # Pathogenic MEI index
    pathogenic_mei_index = defaultdict(list)

    for x in pathogenic_mei.iter_rows(named=True):

        pathogenic_mei_index[x["CHROM"]].append(
            {
                "pos": x["POS"],
                "alt": x[pathogenic_te],
                "annotation":
                    f'Group:"{x[pathogenic_group]}"|'
                    f'Position:"{x[pathogenic_pos]}"|'
                    f'Gene:"{x[pathogenic_gene]}"|'
                    f'Transcript:"{x[pathogenic_transcript]}"|'
                    f'Genomic Feature:"{x[pathogenic_feature]}"|'
                    f'Insertion:"{x[pathogenic_te]}"|'
                    f'Subfamily:"{x[pathogenic_subfamily]}"|'
                    f'Strand:"{x[pathogenic_strand]}"|'
                    f'Size:"{x[pathogenic_size]}"|'
                    f'Disease:"{x[pathogenic_disease]}"|'
                    f'MOI:"{x[pathogenic_moi]}"|'
                    f'PMID:"{x[pathogenic_pmid]}"'
            }
        )

    # ---------------------------------------------------------
    # ClinGen
    # ---------------------------------------------------------
    LOGGER.debug("Loading ClinGen data: %s", args.clingen)

    clingen = (
        pl.read_csv(
            args.clingen,
            skip_rows=4,
            columns=["GENE SYMBOL", "DISEASE LABEL", "MOI", "CLASSIFICATION"],
            schema_overrides={
                "GENE SYMBOL"    : pl.String,
                "DISEASE LABEL"  : pl.String,
                "MOI"            : pl.String,
                "CLASSIFICATION" : pl.String,
            },
        )
        .slice(1)
        .with_columns(
            [
                pl.col("GENE SYMBOL").str.strip_chars(),
                pl.col("DISEASE LABEL").fill_null("").str.strip_chars(),
                pl.col("MOI").fill_null("").str.strip_chars(),
                pl.col("CLASSIFICATION").fill_null("").str.strip_chars(),
            ]
        )
        # Keep only requested ClinGen classifications
        .filter(
            pl.col("CLASSIFICATION").is_in(
                args.clingen_classification
            )
        )
        # Build disease + MOI annotation
        .with_columns(
            (
                pl.col("DISEASE LABEL")
                + " ("
                + pl.col("MOI")
                + ")"
            ).alias("clingen annotation")
        )
        .group_by("GENE SYMBOL")
        .agg(
            pl.col("clingen annotation")
            .unique()
            .sort()
            .str.join("|")
            .alias("clingen annotation")
        )
    )

    clingen_dict = dict(
        zip(
            clingen["GENE SYMBOL"].to_list(),
            clingen["clingen annotation"].to_list(),
        )
    )

    # ---------------------------------------------------------
    # gnomad pLI
    # ---------------------------------------------------------
    LOGGER.debug("Loading gnomAD constraint metrics: %s", args.gnomad_constraint)

    gnomad_constraint = (
        pl.read_csv(
            args.gnomad_constraint,
            separator="\t",
            columns=["gene", "transcript", "lof.pLI"],
            schema_overrides={
                "gene": pl.String,
                "transcript": pl.String,
                "lof.pLI": pl.String,
            },
        )
        .with_columns(
            pl.col("transcript")
            .str.split(".")
            .list.first()
            .alias("transcript")
        )
    )

    gnomad_constraint_dict = dict(
        zip(
            gnomad_constraint["transcript"].to_list(),
            gnomad_constraint["lof.pLI"].to_list(),
        )
    )

    # ---------------------------------------------------------
    # Control cohort
    # ---------------------------------------------------------
    LOGGER.debug("Loading gnomAD control VCFs from %s", args.control_dir)

    control_af_index = defaultdict(list)

    for control_vcf in control_vcfs:

        vcf = open_vcf(
            os.path.join(
                args.control_dir,
                control_vcf,
            )
        )

        with redirect_htslib_stderr():
            for v in vcf:

                af = v.INFO.get(control_af)

                pop_af = {
                    pop: v.INFO.get(field)
                    for field, pop in POP_AF_FIELDS.items()
                }
                pop_af = {p: float(a) for p, a in pop_af.items() if a is not None}

                if pop_af:
                    popmax_pop, popmax_af = max(pop_af.items(), key=lambda x: x[1])
                else:
                    popmax_pop, popmax_af = ".", None

                for alt in v.ALT:
                    control_af_index[v.CHROM].append({
                        "pos": v.POS,
                        "alt": alt,
                        "af": af,
                        "grpmax_af": popmax_af,
                        "grpmax_pop": popmax_pop,
                    })

    # ---------------------------------------------------------
    # Process VCFs
    # ---------------------------------------------------------

    all_files = []
    qc_failed_by_file = {}

    for original_vcf_file in vcf_files:

        LOGGER.debug("Processing VCF: %s", original_vcf_file)

        prepared_vcf = prepare_vcf(
            original_vcf_file,
            harmonize_te
        )

        vcf = open_vcf(prepared_vcf)

        LOGGER.debug("Running MANE intersection")
        with redirect_htslib_stderr():
            mane_intersection_dict = mane_intersection(
                args.mane,
                vcf,
                mode,
                args.alt
            )

        # Need a fresh VCF because the first iteration consumed it
        vcf = open_vcf(prepared_vcf)

        with redirect_htslib_stderr():
            annotations, qc_failed_count = annotate_vcf(
                vcf,
                mane_intersection_dict,
                args.alt,
                source=os.path.basename(original_vcf_file),
                clingen_only=clingen_only,
                qc=args.qc,
                siteonly=siteonly,
                window=window,
                pathogenic_mei_index=pathogenic_mei_index,
                control_af_index=control_af_index,
                clingen_dict=clingen_dict,
                gnomad_constraint_dict=gnomad_constraint_dict,
            )

        all_files.extend(annotations)

        qc_failed_by_file[
            os.path.basename(original_vcf_file)
        ] = qc_failed_count

    # ---------------------------------------------------------
    # Export
    # order
    # 1. pathogenic MEI
    # 2. CDS > junction > UTR5 > UTR3
    # 3. lof.pLI descending
    # ---------------------------------------------------------

    if not all_files:
        LOGGER.debug("No variants passed QC/filtering.")
        return qc_failed_by_file

    result = pl.DataFrame(all_files)
    
    result = result.with_columns(
        [
            pl.col("mapped_pathogenic_MEI")
            .eq(".")
            .cast(pl.Int8)
            .alias("_pathogenic_MEI_sort"),

            pl.col("feature")
            .map_elements(
                feature_rank,
                return_dtype=pl.Int64,
            )
            .alias("_feature_sort"),

            pl.col("lof.pLI")
            .cast(pl.Float64, strict=False)
            .alias("_lof_pLI_sort"),
        ]
    )

    result = result.sort(
        [
            "_pathogenic_MEI_sort",
            "_feature_sort",
            "_lof_pLI_sort",
        ],
        descending=[
            False,
            False,
            True,
        ],
        nulls_last=True,
    )

    result = (
        result
        .drop([
            "_pathogenic_MEI_sort",
            "_feature_sort",
            "_lof_pLI_sort",
        ])
        .with_columns(
            pl.all()
            .cast(pl.String)
            .fill_null(".")
            .replace("", ".")
        )
    )

    result.write_csv(
        args.output,
        include_header=True,
    )

    LOGGER.debug("Wrote %d rows to %s", result.height, args.output)

    return qc_failed_by_file


if __name__ == "__main__":
    main()