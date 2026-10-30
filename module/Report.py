import argparse
import base64
import html
import re
import math
import logging
from pathlib import Path
import polars as pl

LOGGER = logging.getLogger(__name__)

# ============================================================
# Configuration
# ============================================================

KNOWN_COLUMNS = {
    "CHROM",
    "POS",
    "ME",
    "ALT",
    "gene",
    "feature",
    "SVLEN",
    "strand",
    "mapped_pathogenic_MEI",
    "gnomad_AF",
    "clingen annotation",
    "lof.pLI",
    "sample",
}

INTERNAL_COLUMNS = {
    "_af_num",
    "_pli_num",
    "_af_zero",
    "_dominant",
    "_moi",
    "_moi_rank",
}

# ============================================================
# Utility functions
# ============================================================

def safe(v):
    """Convert None/NaN to a printable string."""
    if v is None:
        return "."

    if isinstance(v, float) and math.isnan(v):
        return "."

    return str(v)


def esc(v):
    return html.escape(safe(v))


def number(v, default=None):
    try:
        if v is None:
            return default

        if isinstance(v, float) and math.isnan(v):
            return default

        if v == "." or v == "":
            return default

        return float(v)

    except (ValueError, TypeError):
        return default


def split_pipe(s):
    """Split 'A|B|C' while dropping empty values."""
    return [x.strip() for x in safe(s).split("|") if x.strip() and x.strip() != "."]


def parse_kv_string(s):
    """
    Parse:
        Group:"Mendelian disorders"|Position:"chr4:15502903"|Strand:"+"|...

    into:
        {
            "Group": "Mendelian disorders",
            "Position": "chr4:15502903",
            ...
        }
    """
    result = {}

    for item in split_pipe(s):
        if ":" not in item:
            continue

        key, value = item.split(":", 1)
        key = key.strip()
        value = value.strip()

        # Remove surrounding quotes.
        if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
            value = value[1:-1]

        result[key] = value

    return result

def parse_gnomad_af(v):
    """
    Parse:
        chrx:xxx|<INS:ME:xxx>|xxx|MAF:xxx|MAF_POP:xx
    """

    parts = [x.strip() for x in safe(v).split("|")]

    result = {
        "position": parts[0] if len(parts) > 0 else ".",
        "mei": parts[1].strip("<>").split(":")[-1] if len(parts) > 1 else ".",
        "af": number(parts[2]) if len(parts) > 2 else None,
        "maf": None,
        "maf_pop": ".",
    }

    for part in parts[3:]:
        if ":" not in part:
            continue

        key, value = part.split(":", 1)

        if key == "MAF":
            result["maf"] = number(value)
        elif key == "MAF_POP":
            result["maf_pop"] = value


    return result


def parse_moi_from_pathogenic(v):
    kv = parse_kv_string(v)
    return kv.get("MOI", ".")


def is_common_polymorphism(v):
    group = parse_kv_string(v).get("Group", ".")
    return group == "Common polymorphism"


def clingen_has_dominant(annotation):
    """
    AD/XLD gets highest priority.
    """
    s = safe(annotation)

    # Match standalone AD or XLD.
    return bool(re.search(r"\b(?:AD|XLD)\b", s, re.I))


def format_pli(v):
    x = number(v)
    if x is None:
        return "."
    return "0" if x == 0 else f"{x:.3g}"


def format_af(af):

    if af is None:
        return "."

    if af == 0: # prevent 0.0
        return "0"

    return f"{af:.3g}"

def format_feature(feature):

    return feature.replace(".", "intron").replace("_", " ")


def pli_color_score(v):
    """
    pLI:
        >0.8 -> red
        low   -> green
    """
    x = number(v)

    if x is None:
        return 0.0

    return max(0.0, min(1.0, (x - 0.0) / 1.0))


def af_color_score(v):
    """
    AF:
        missing/0      -> strongest red
        extremely rare -> red
        increasing AF  -> increasingly green
    """
    x = number(v)

    if x is None or x == 0:
        return 1.0

    import math

    # Log-scale because AF spans many orders of magnitude.
    return max(
        0.0,
        min(
            1.0,
            1.0 - (math.log10(x) + 6) / 6,
        ),
    )


def color_style(score):
    """
    Morandi-style red -> green.

    score = 1 => red
    score = 0 => green
    """
    score = max(0.0, min(1.0, score))

    # Muted/Morandi colors.
    red = (205, 150, 145)
    green = (160, 185, 165)

    r = int(green[0] * (1 - score) + red[0] * score)
    g = int(green[1] * (1 - score) + red[1] * score)
    b = int(green[2] * (1 - score) + red[2] * score)

    return f"background: rgb({r},{g},{b});"


def moi_rank(moi):
    """
    Red -> green priority:
        AD > XLD > AR > XLR > all others
    """
    s = safe(moi)

    if re.search(r"\bAD\b", s, re.I):
        return 0

    if re.search(r"\bXLD\b", s, re.I):
        return 1

    if re.search(r"\bAR\b", s, re.I):
        return 2

    if re.search(r"\bXLR\b", s, re.I):
        return 3

    return 4


def find_plot(plot_dir, chrom, pos, alt):
    """
    Expected:
        {CHROM}_{POS}_{ALT}.png
    """
    filename = f"{chrom}_{pos}_{alt}.png"
    path = Path(plot_dir) / filename

    if path.exists():
        return path

    return None


def png_as_data_uri(path):
    if path is None:
        return None

    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"

def sample_gt_html(value):
    """
    Convert the sample column into readable sample:GT entries.

    Example input:
        sample1:0/1;sample2:0/0;sample3:1/1

    Output:
        sample1:0/1; sample2:0/0; sample3:1/1
    """
    if value is None or (
        isinstance(value, float) and math.isnan(value)
    ) or str(value).strip() in {"", "."}:
        return "."

    return "; ".join(
        x.strip()
        for x in str(value).split(";")
        if x.strip()
    )

# ============================================================
# HTML rendering
# ============================================================

HTML_TEMPLATE = r"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>MEander Report</title>

<style>
    :root {
        --bg: #f4f2ef;
        --card: #ffffff;
        --text: #383838;
        --muted: #77736f;
        --border: #d8d3ce;
        --accent: #827b75;
        --common: #edf5ee;
        --common-border: #cbdcca;
    }

    * {
        box-sizing: border-box;
    }

    body {
        margin: 0;
        padding: 0;
        background: var(--bg);
        color: var(--text);
        font-family:
            -apple-system,
            BlinkMacSystemFont,
            "Segoe UI",
            Arial,
            sans-serif;
    }

    .header {
        padding: 28px 36px 18px;
        background: #e8e4df;
        border-bottom: 1px solid var(--border);
    }

    h1 {
        margin: 0 0 8px;
        font-size: 26px;
    }

    .subtitle {
        color: var(--muted);
        font-size: 14px;
    }

    .tabs {
        display: flex;
        gap: 8px;
        padding: 18px 36px 0;
    }

    .tab-button {
        border: 1px solid var(--border);
        background: #ddd9d4;
        color: var(--text);
        padding: 10px 16px;
        border-radius: 8px 8px 0 0;
        cursor: pointer;
        font-size: 14px;
    }

    .tab-button.active {
        background: var(--card);
        border-bottom-color: var(--card);
    }

    .page {
        display: none;
        padding: 20px 36px 40px;
    }

    .page.active {
        display: block;
    }

    details.mei {
        background: var(--card);
        border: 1px solid var(--border);
        border-radius: 10px;
        margin: 10px 0;
        overflow: hidden;
    }

    details.mei.common-polymorphism {
        background: var(--common);
        border-color: var(--common-border);
    }

    details.mei > summary {
        cursor: pointer;
        list-style: none;
        padding: 13px 16px;
    }

    details.mei.common-polymorphism > summary {
        background: var(--common);
    }

    details.mei > summary::-webkit-details-marker {
        display: none;
    }

    details.mei > summary::before {
        content: "▶";
        display: inline-block;
        margin-right: 10px;
        font-size: 11px;
        color: var(--muted);
    }

    details.mei[open] > summary::before {
        content: "▼";
    }

    .summary-grid {
        display: grid;
        grid-template-columns:
            minmax(120px, 0.7fr)
            minmax(100px, 0.6fr)
            minmax(120px, 0.9fr)
            minmax(180px, 1.4fr)
            minmax(120px, 0.9fr)
            minmax(80px, 0.6fr)
            minmax(100px, 0.7fr);
        gap: 8px;
        align-items: stretch;
    }

    .cell {
        padding: 8px 10px;
        border-radius: 6px;
        min-width: 0;
    }

    .small-cell {
        padding: 8px 8px;
        border-radius: 6px;
        min-width: 0;
    }

    .label {
        font-size: 10px;
        text-transform: uppercase;
        letter-spacing: .05em;
        color: var(--muted);
        margin-bottom: 3px;
    }

    .value {
        font-size: 13px;
        overflow-wrap: anywhere;
    }

    .content {
        border-top: 1px solid var(--border);
        padding: 18px;
    }

    .content-grid {
        display: grid;
        grid-template-columns: minmax(300px, 1fr) minmax(300px, 1fr);
        gap: 20px;
    }

    .plot {
        width: 100%;
        max-height: 700px;
        object-fit: contain;
        border: 1px solid var(--border);
        border-radius: 6px;
        background: #faf9f7;
    }

    .kv-table {
        width: 100%;
        border-collapse: collapse;
        font-size: 12px;
    }

    .kv-table th,
    .kv-table td {
        text-align: left;
        padding: 7px 8px;
        border-bottom: 1px solid #e8e4df;
        vertical-align: top;
        overflow-wrap: anywhere;
    }

    .kv-table th {
        width: 180px;
        color: var(--muted);
        font-weight: 500;
    }

    .annotation {
        display: inline-block;
        min-width: 50px;
        margin: 2px 4px 2px 0;
        padding: 4px 7px;
        border-radius: 5px;
        background: #ece9e5;
        font-size: 12px;
        text-align: center;
    }

    .moi {
        font-weight: 600;
        padding: 4px 8px;
        border-radius: 5px;
        display: inline-block;
    }

    .empty {
        color: var(--muted);
        padding: 20px;
        text-align: center;
    }

    @media (max-width: 1000px) {
        .summary-grid {
            grid-template-columns: repeat(2, 1fr);
        }

        .content-grid {
            grid-template-columns: 1fr;
        }

        .page,
        .tabs,
        .header {
            padding-left: 15px;
            padding-right: 15px;
        }
    }
</style>

<script>
function showPage(pageId, button) {
    document.querySelectorAll(".page").forEach(
        x => x.classList.remove("active")
    );

    document.querySelectorAll(".tab-button").forEach(
        x => x.classList.remove("active")
    );

    document.getElementById(pageId).classList.add("active");
    button.classList.add("active");
}
</script>
</head>

<body>

<div class="header">
    <h1>MEander Report</h1>
    <div class="subtitle">
        {{ gene_count }} MEIs mapped to coding genes ·
        {{ pathogenic_count }} MEIs mapped to pathogenic MEIs
    </div>
</div>

<div class="tabs">
    <button
        class="tab-button active"
        onclick="showPage('gene-page', this)">
        MEI → Coding genes ({{ gene_count }})
    </button>

    <button
        class="tab-button"
        onclick="showPage('pathogenic-page', this)">
        MEI → Pathogenic MEI ({{ pathogenic_count }})
    </button>
</div>

<div id="gene-page" class="page active">
    {{ gene_sections }}
</div>

<div id="pathogenic-page" class="page">
    {{ pathogenic_sections }}
</div>

</body>
</html>
"""


# ============================================================
# Section rendering
# ============================================================

def render_annotations(annotation):
    annotations = split_pipe(annotation)

    if not annotations:
        return "."

    output = []

    for ann in annotations:

        # AD/XLD annotations get red-side highlighting.
        rank = moi_rank(ann)
        style = color_style(1.0 - (rank / 4.0)) if rank < 4 else ""

        output.append(
            f'<span class="annotation" style="{style}">{esc(ann)}</span>'
        )

    return "".join(output)


def render_plot(plot_path):
    if plot_path is None:
        return '<div class="empty">Plot not found</div>'

    uri = png_as_data_uri(plot_path)

    return f'<img class="plot" src="{uri}" alt="MEI plot">'


def render_additional_columns(row, excluded_columns):
    rows = []

    for col, value in row.items():

        if col in excluded_columns:
            continue
        if col in INTERNAL_COLUMNS:
            continue

        value = row[col]

        if value is None:
            continue
        if isinstance(value, float) and math.isnan(value):
            continue

        rows.append(
            f"""
            <tr>
                <th>{esc(col)}</th>
                <td>{esc(value)}</td>
            </tr>
            """
        )

    if not rows:
        return ""

    return f"""
    <table class="kv-table">
        {''.join(rows)}
    </table>
    """

def render_gnomad_af(v):
    parsed = parse_gnomad_af(v)

    return f"""
    <div>
        <i>POS</i>: {esc(parsed["position"])}
        &nbsp;&nbsp;
        <i>ALT</i>: {esc(parsed["mei"])}
        &nbsp;&nbsp;
        <i>AF</i>: {esc(format_af(parsed["af"]))}
        &nbsp;&nbsp;
        <i>MAF</i>: {esc(format_af(parsed["maf"]))}
        &nbsp;&nbsp;
        <i>MAF population</i>: {esc(parsed["maf_pop"])}
    </div>
    """


def render_gene_section(row, plot_dir):
    chrom = safe(row.get("CHROM"))
    pos = safe(row.get("POS"))
    alt = safe(row.get("ALT"))

    maf = parse_gnomad_af(row.get("gnomad_AF"))["maf"]
    pli = number(row.get("lof.pLI"))

    # gnomAD: zero => red.
    af_score_ = af_color_score(maf)
    af_style = color_style(af_score_)

    # pLI > 0.8 => red, otherwise gradient.
    pli_score_ = pli_color_score(pli)
    pli_style = color_style(pli_score_)

    clingen = safe(row.get("clingen annotation"))
    feature = safe(format_feature(row.get("feature")))

    plot = find_plot(plot_dir, chrom, pos, alt)

    excluded = {
        "CHROM",
        "POS",
        "ME",
        "ALT",
        "gene",
        "clingen annotation",
        "lof.pLI",
        "gnomad_AF",
        "mapped_pathogenic_MEI",
        "sample",
    }

    additional = render_additional_columns(row, excluded)

    return f"""
    <details class="mei">
        <summary>
            <div class="summary-grid">

                <div class="small-cell">
                    <div class="label">Breakpoint</div>
                    <div class="value">{esc(chrom)}:{esc(pos)}</div>
                </div>

                <div class="small-cell">
                    <div class="label">MEI</div>
                    <div class="value">{esc(row.get("ALT", alt))}</div>
                </div>

                <div class="cell">
                    <div class="label">Gene</div>
                    <div class="value">{esc(row.get("gene"))}</div>
                </div>

                <div class="cell">
                    <div class="label">ClinGen annotation</div>
                    <div class="value">{render_annotations(clingen)}</div>
                </div>

                <div class="cell">
                    <div class="label">Feature</div>
                    <div class="value">{esc(feature)}</div>
                </div>

                <div class="cell">
                    <div class="label">pLI</div>
                    <div class="value">
                        <span class="annotation" style="{pli_style}">
                            {esc(format_pli(row.get("lof.pLI")))}
                        </span>
                    </div>
                </div>

                <div class="cell">
                    <div class="label">gnomAD MAF</div>
                    <div class="value">
                        <span class="annotation" style="{af_style}">
                            {esc(format_af(maf))}
                        </span>
                    </div>
                </div>

            </div>
        </summary>

        <div class="content">
            <div class="content-grid">

                <div>
                    {render_plot(plot)}
                </div>

                <div>
                    <table class="kv-table">

                        <tr>
                            <th>Breakpoint position</th>
                            <td>{esc(chrom)}:{esc(pos)}</td>
                        </tr>

                        <tr>
                            <th>MEI</th>
                            <td>{esc(row.get("ALT", alt))}</td>
                        </tr>

                        <tr>
                            <th>Gene</th>
                            <td>{esc(row.get("gene"))}</td>
                        </tr>

                        <tr>
                            <th>ClinGen annotation</th>
                            <td>{render_annotations(clingen)}</td>
                        </tr>

                        <tr>
                            <th>lof.pLI</th>
                            <td>{esc(format_pli(row.get("lof.pLI")))}</td>
                        </tr>

                        <tr>
                            <th>gnomAD</th>
                            <td>{render_gnomad_af(row.get("gnomad_AF"))}</td>
                        </tr>

                        <tr>
                            <th>Samples</th>
                            <td>{esc(sample_gt_html(row.get("sample")))}</td>
                        </tr>

                    </table>

                    {additional}
                </div>

            </div>
        </div>
    </details>
    """


def render_pathogenic_section(row, plot_dir):
    chrom = safe(row.get("CHROM"))
    pos = safe(row.get("POS"))
    alt = safe(row.get("ALT"))

    mapped = parse_kv_string(row.get("mapped_pathogenic_MEI"))
    moi = mapped.get("MOI", ".")

    # AD red -> XLD -> AR -> XLR -> others green.
    rank = moi_rank(moi)
    moi_score = 1.0 - (rank / 4.0)
    moi_style = color_style(moi_score)

    maf = parse_gnomad_af(row.get("gnomad_AF"))["maf"]
    pli = number(row.get("lof.pLI"))

    af_score_ = af_color_score(maf)
    af_style = color_style(af_score_)

    pli_score_ = pli_color_score(pli)
    pli_style = color_style(pli_score_)

    clingen = safe(row.get("clingen annotation"))
    feature = safe(format_feature(row.get("feature")))

    plot = find_plot(plot_dir, chrom, pos, alt)

    excluded = {
        "CHROM",
        "POS",
        "ME",
        "ALT",
        "gene",
        "clingen annotation",
        "lof.pLI",
        "gnomad_AF",
        "mapped_pathogenic_MEI",
        "sample",
    }

    additional = render_additional_columns(row, excluded)

    mapped_rows = []

    for key, value in mapped.items():
        mapped_rows.append(
            f"""
            <tr>
                <th>{esc(key)}</th>
                <td>{esc(value)}</td>
            </tr>
            """
        )

    mapped_html = "".join(mapped_rows)

    group = mapped.get("Group", ".")
    card_class = (
        "mei common-polymorphism"
        if group == "Common polymorphism"
        else "mei"
    )

    return f"""
    <details class="{card_class}">
        <summary>
            <div class="summary-grid">

                <div class="cell">
                    <div class="label">Breakpoint</div>
                    <div class="value">{esc(chrom)}:{esc(pos)}</div>
                </div>

                <div class="cell">
                    <div class="label">MEI</div>
                    <div class="value">{esc(row.get("ALT", alt))}</div>
                </div>

                <div class="cell">
                    <div class="label">Gene</div>
                    <div class="value">{esc(row.get("gene"))}</div>
                </div>

                <div class="cell">
                    <div class="label">Feature</div>
                    <div class="value">{esc(feature)}</div>
                </div>

                <div class="cell">
                    <div class="label">Mode of inheritance</div>
                    <div class="value">
                        <span class="annotation" style="{moi_style}">
                            {esc(moi)}
                        </span>
                    </div>
                </div>

                <div class="cell">
                    <div class="label">pLI</div>
                    <div class="value">
                        <span class="annotation" style="{pli_style}">
                            {esc(format_pli(row.get("lof.pLI")))}
                        </span>
                    </div>
                </div>

                <div class="cell">
                    <div class="label">gnomAD MAF</div>
                        <span class="annotation" style="{af_style}">
                            {esc(format_af(maf))}
                        </span>
                </div>

            </div>
        </summary>

        <div class="content">
            <div class="content-grid">

                <div>
                    {render_plot(plot)}
                </div>

                <div>
                    <table class="kv-table">

                        <tr>
                            <th>Breakpoint position</th>
                            <td>{esc(chrom)}:{esc(pos)}</td>
                        </tr>

                        <tr>
                            <th>MEI</th>
                            <td>{esc(row.get("ALT", alt))}</td>
                        </tr>

                        <tr>
                            <th>Gene</th>
                            <td>{esc(row.get("gene"))}</td>
                        </tr>

                        <tr>
                            <th>Mode of inheritance</th>
                            <td>
                                <span class="moi" style="{moi_style}">
                                    {esc(moi)}
                                </span>
                            </td>
                        </tr>

                        <tr>
                            <th>Mapped disease-relevant MEI</th>
                            <td>
                                <table class="kv-table">
                                    {mapped_html}
                                </table>
                            </td>
                        </tr>

                        <tr>
                            <th>ClinGen annotation</th>
                            <td>
                                {render_annotations(clingen)}
                            </td>
                        </tr>

                        <tr>
                            <th>lof.pLI</th>
                            <td>{esc(format_pli(row.get("lof.pLI")))}</td>
                        </tr>

                        <tr>
                            <th>gnomAD</th>
                            <td>{render_gnomad_af(row.get("gnomad_AF"))}</td>
                        </tr>

                        <tr>
                            <th>Samples</th>
                            <td>{esc(sample_gt_html(row.get("sample")))}</td>
                        </tr>

                    </table>

                    {additional}
                </div>

            </div>
        </div>
    </details>
    """


# ============================================================
# Main
# ============================================================

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate an interactive MEI HTML report."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input TSV/CSV table.",
    )

    parser.add_argument(
        "-p",
        "--plots",
        required=True,
        help="Directory containing {CHROM}_{POS}_{ME}.png files.",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="mei_report.html",
        help="Output HTML file.",
    )

    parser.add_argument(
        "--no-common",
        action="store_true",
        help="Exclude MEIs mapped to common polymorphism group.",
    )

    args = parser.parse_args(argv)

    # --------------------------------------------------------
    # Read table
    # --------------------------------------------------------

    input_path = Path(args.input)

    suffix = input_path.suffix.lower()

    if suffix == ".tsv":
        sep = "\t"
    elif suffix == ".csv":
        sep = ","
    else:
        raise ValueError(
            f"Unsupported input format: {input_path.suffix!r}. "
            "Expected .tsv or .csv."
        )

    df = pl.read_csv(
        args.input,
        separator=sep,
        infer_schema=False,
        null_values=[],
    )

    required = [
        "CHROM",
        "POS",
        "ALT",
        "mapped_pathogenic_MEI",
        "gnomad_AF",
        "clingen annotation",
        "lof.pLI",
    ]

    missing = [x for x in required if x not in df.columns]

    if missing:
        raise ValueError(
            "Missing required columns: " + ", ".join(missing)
        )

    df = df.with_columns(
        pl.col("gnomad_AF")
        .map_elements(
            lambda x: parse_gnomad_af(x)["af"],
            return_dtype=pl.Float64,
        )
        .alias("_af_num"),

        pl.col("lof.pLI")
        .cast(pl.Float64, strict=False)
        .alias("_pli_num"),
    )

    # --------------------------------------------------------
    # Separate MEI types
    # --------------------------------------------------------

    # if a MEI mapped to both pathogenic and gene, assign to pathogenic group
    pathogenic_mask = (
        pl.col("mapped_pathogenic_MEI")
        .cast(pl.String)
        .str.strip_chars()
        .ne("")
        & pl.col("mapped_pathogenic_MEI")
        .cast(pl.String)
        .str.strip_chars()
        .ne(".")
    )

    pathogenic_df = df.filter(pathogenic_mask)
    gene_df = df.filter(~pathogenic_mask)

    # --------------------------------------------------------
    # Optionally exclude common polymorphisms
    # --------------------------------------------------------

    if args.no_common:
        pathogenic_df = pathogenic_df.filter(
            ~pl.col("mapped_pathogenic_MEI").map_elements(
                is_common_polymorphism,
                return_dtype=pl.Boolean,
            )
        )

    # --------------------------------------------------------
    # Sorting
    # --------------------------------------------------------

    # Gene page:
    #   1. gnomAD AF == 0
    #   2. high pLI
    #   3. ClinGen AD/XLD
    #
    # These are deliberately sorted in the same direction as
    # the visual highlighting: strongest candidates first.

    if gene_df.height > 0:
        gene_df = gene_df.with_columns(
            pl.col("_af_num")
            .eq(0)
            .fill_null(False)
            .alias("_af_zero"),

            pl.col("clingen annotation")
            .map_elements(
                clingen_has_dominant,
                return_dtype=pl.Boolean,
            )
            .alias("_dominant"),
        )

        gene_df = gene_df.sort(
            by=[
                "_af_zero",
                "_pli_num",
                "_dominant",
            ],
            descending=[
                True,
                True,
                True,
            ],
            nulls_last=True,
        )

    # Pathogenic page:
    # AD > XLD > AR > XLR > others

    if pathogenic_df.height > 0:
        pathogenic_df = pathogenic_df.with_columns(
            pl.col("mapped_pathogenic_MEI")
            .map_elements(
                parse_moi_from_pathogenic,
                return_dtype=pl.String,
            )
            .alias("_moi"),
        )

        pathogenic_df = pathogenic_df.with_columns(
            pl.col("_moi")
            .map_elements(
                moi_rank,
                return_dtype=pl.Int64,
            )
            .alias("_moi_rank"),
        )

        pathogenic_df = pathogenic_df.sort(
            by=[
                "_moi_rank",
                "_af_num",
                "_pli_num",
            ],
            descending=[
                False,
                False,
                True,
            ],
            nulls_last=True,
        )
    # --------------------------------------------------------
    # Render sections
    # --------------------------------------------------------

    gene_sections = []

    for row in gene_df.iter_rows(named=True):
        gene_sections.append(
            render_gene_section(
                row,
                args.plots,
            )
        )

    pathogenic_sections = []

    for row in pathogenic_df.iter_rows(named=True):
        pathogenic_sections.append(
            render_pathogenic_section(
                row,
                args.plots,
            )
        )

    gene_html = (
        "\n".join(gene_sections)
        if gene_sections
        else '<div class="empty">No MEIs mapped to coding genes.</div>'
    )

    pathogenic_html = (
        "\n".join(pathogenic_sections)
        if pathogenic_sections
        else '<div class="empty">No MEIs mapped to pathogenic MEIs.</div>'
    )

    # --------------------------------------------------------
    # Fill template
    # --------------------------------------------------------

    html_output = (
        HTML_TEMPLATE
        .replace("{{ gene_count }}", str(gene_df.height))
        .replace("{{ pathogenic_count }}", str(pathogenic_df.height))
        .replace("{{ gene_sections }}", gene_html)
        .replace("{{ pathogenic_sections }}", pathogenic_html)
    )

    Path(args.output).write_text(
        html_output,
        encoding="utf-8",
    )

    print(f"Generated: {args.output}")
    print(f"Gene-feature MEIs: {gene_df.height}")
    print(f"Pathogenic-MEI mappings: {pathogenic_df.height}")


if __name__ == "__main__":
    main()
