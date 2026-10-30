from __future__ import annotations
__version__ = "0.1.0"

import argparse
import logging
import warnings
import sys
from pathlib import Path

from module import Mapping
from module import Drawing
from module import Report

"""
MEander entry.

Mapping_module.py -> [annotated table] -> Drawing_module.py -> [plots] -> Report_module.py -> [HTML report]
"""

LOGGER = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def configure_logging(log_file: Path, verbose: bool = False) -> None:
    """Configure application-wide logging."""
    level = logging.DEBUG if verbose else logging.INFO

    logging.basicConfig(
        filename=log_file,
        level=level,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logging.getLogger("PIL").setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)


# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        prog="MEander",
        description="Run the complete MEI mapping -> drawing -> reporting pipeline."
    )

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Print version information",
    )

    # General
    general = parser.add_argument_group("general options")

    general.add_argument(
        "--input",
        required=True,
        type=Path,
        help="[REQUIRED] VCF or list of VCFs. If list, the file must end with .list.txt with one VCF path per line",
    )

    general.add_argument(
        "--output",
        required=True,
        help="[REQUIRED] Output prefix",
    )

    general.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="[REQUIRED] Directory for MEander outputs",
    )

    general.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose/debug logging",
    )

    # Mapping
    mapping = parser.add_argument_group("mapping options")

    mapping.add_argument(
        "--position-to-map",
        choices=[
                "start",
                "end",
                ],
        default="start",
        help=(
                "Position in the disease-associated MEI catalog "
                "(TSD start or end) to use. For MELT input, "
                "end is recommended (default: start)"
            ),
    )

    mapping.add_argument(
        "--mane",
        required=True,
        type=Path,
        help="[REQUIRED] MANE annotation BED file",
    )

    mapping.add_argument(
        "--pathogenic",
        required=True,
        type=Path,
        help="[REQUIRED] Disease-focused MEI catalog file",
    )

    mapping.add_argument(
        "--clingen",
        required=True,
        type=Path,
        help="[REQUIRED] ClinGen gene-disease summary file",
    )

    mapping.add_argument(
        "--gnomad-constraint",
        required=True,
        type=Path,
        help="[REQUIRED] gnomAD constraint metrics file",
    )

    mapping.add_argument(
        "--control-dir",
        required=True,
        type=Path,
        help="[REQUIRED] Directory with gnomAD resources",
    )

    mapping.add_argument(
        "--mode",
        nargs="+",
        choices=[
            "CDS",
            "junction_donor",
            "junction_acceptor",
            "three_prime_UTR",
            "five_prime_UTR",
        ],
        metavar="MODE",
        default=["CDS"],
        help=(
            "Keep MEI in selected genomic features.  "
            "[CDS, junction_donor, junction_acceptor, "
            "three_prime_UTR, five_prime_UTR]"
        ),
    )

    mapping.add_argument(
        "--clingen-only",
        action="store_true",
        help="Keep MEI in genes with ClinGen Gene-Disease Clinical Validity",
    )

    mapping.add_argument(
        "--qc",
        choices=[
            "MELT",
            "MEGAnE",
            "xTea",
        ],
        default=None,
        help=(
            "Remove low-confidence call. "
            "MELT: ASSESS >=3, split reads >2, FILTER=PASS; "
            "MEGAnE: FILTER=PASS. "
            "(default: none)"
        ),
    )

    mapping.add_argument(
        "--siteonly",
        action="store_true",
        help="When enabled, sample-level information is not collected",
    )

    mapping.add_argument(
        "--window",
        type=int,
        default=20,
        help="Mapping window size (default: 20)",
    )

    mapping.add_argument(
        "--alt",
        default="ALT",
        help="Field with MEI family (default: ALT)",
    )

    mapping.add_argument(
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
        metavar="CLASSIFICATION",
        default=["Definitive", "Strong"],
        help=(
            "ClinGen classifications to include.  "
            "[Definitive, Strong, "
            "Moderate, Limited, No Known Disease Relationship, "
            "Disputed, Refuted]"
        ),
    )

    # Drawing
    drawing = parser.add_argument_group("drawing options")

    drawing.add_argument(
        "--gff",
        required=True,
        type=Path,
        help="[REQUIRED] GFF3 annotation file",
    )

    drawing.add_argument(
        "--melt-meinfo",
        action="store_true",
        help="Handle MEI length and strand in MELT MEINFO explicitly",
    )

    drawing.add_argument(
        "--melen-field",
        default=None,
        help="VCF field or INFO tag with MEI length",
    )

    drawing.add_argument(
        "--mestrand-field",
        default=None,
        help="VCF field or INFO tag with MEI strand",
    )

    # Report
    report = parser.add_argument_group("report options")

    report.add_argument(
        "--no-common",
        action="store_true",
        help="Exclude common variants from the report",
    )

    return parser


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def run_mapping(args: argparse.Namespace, annotated_csv: Path) -> None:
    """Run Mapping_module."""
    LOGGER.info("Step 1/3: running mapping")        

    mapping_args = [
        "--input", str(args.input),
        "--output", str(annotated_csv),
        "--mane", str(args.mane),
        "--pathogenic", str(args.pathogenic),
        "--clingen", str(args.clingen),
        "--gnomad-constraint", str(args.gnomad_constraint),
        "--control-dir", str(args.control_dir),
        "--mode", *args.mode,
        "--window", str(args.window),
        "--position-to-map", str(args.position_to_map),
        "--alt", args.alt,
        "--clingen-classification", *args.clingen_classification,
    ]

    if args.qc is not None:
        LOGGER.info("QC method: %s", args.qc)
        mapping_args.extend(["--qc", args.qc])
    else:
        LOGGER.info("QC: disabled")

    if args.clingen_only:
        LOGGER.info("Keeping only MEIs in genes with ClinGen GDV")
        mapping_args.append("--clingen-only")

    if args.siteonly:
        LOGGER.info("Treat as site-only VCF")
        mapping_args.append("--siteonly")

    qc_failed_by_file = Mapping.main(mapping_args)

    if not annotated_csv.exists():
        raise RuntimeError(
            f"Mapping completed but expected output was not created: "
            f"{annotated_csv}"
        )

    for vcf_file, count in qc_failed_by_file.items():
        LOGGER.info(
            "QC: %s -> %d variant(s) failed QC",
            vcf_file,
            count,
        )

    LOGGER.info("Mapping completed: %s", annotated_csv)


def run_drawing(args: argparse.Namespace, annotated_csv: Path) -> None:
    """Run Drawing_module."""
    LOGGER.info("Step 2/3: generating plots")

    plots_path = args.output_dir / "plots"

    drawing_args = [
        "--input", str(annotated_csv),
        "--melen-field", str(args.melen_field),
        "--mestrand-field", str(args.mestrand_field),
        "--gff", str(args.gff),
        "--output-dir", str(plots_path),
    ]

    if args.melt_meinfo:
        drawing_args.append("--melt-meinfo")

    Drawing.main(drawing_args)

    LOGGER.info("Drawing completed: %s", args.output_dir)


def run_report(args: argparse.Namespace, annotated_csv: Path) -> Path:
    """Run Report_module."""
    LOGGER.info("Step 3/3: generating report")

    plot_path = args.output_dir / "plots"
    report_path = args.output_dir / f"{args.output}_report.html"

    report_args = [
        "-i", str(annotated_csv),
        "-p", str(plot_path),
        "-o", str(report_path),
    ]

    if args.no_common:
        report_args.append("--no-common")

    Report.main(report_args)

    if not report_path.exists():
        raise RuntimeError(
            f"Report completed but expected output was not created: "
            f"{report_path}"
        )

    LOGGER.info("Report completed: %s", report_path)

    return report_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if not args.melt_meinfo:
        missing = []

        if args.melen_field is None:
            missing.append("--melen-field")

        if args.mestrand_field is None:
            missing.append("--mestrand-field")

        if missing:
            warnings.warn(
                f"{', '.join(missing)} is not provided."
            )

    log_file = args.output_dir / f"{args.output}_meander.log"
    configure_logging(log_file, args.verbose)

    try:
        # Prepare output directory.
        args.output_dir.mkdir(parents=True, exist_ok=True)

        annotated_csv = args.output_dir / f"{args.output}.annotated.csv"

        LOGGER.info("Starting MEander")
        LOGGER.info("Input: %s", args.input)
        LOGGER.info("Output directory: %s", args.output_dir)

        # Pipeline
        run_mapping(args, annotated_csv)
        run_drawing(args, annotated_csv)
        report_path = run_report(args, annotated_csv)

        LOGGER.info("Pipeline completed successfully")
        LOGGER.info("Final report: %s", report_path)

        return 0

    except KeyboardInterrupt:
        LOGGER.warning("Pipeline interrupted by user")
        return 130

    except Exception:
        LOGGER.exception("Pipeline failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
