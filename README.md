# MEander

MEander is a tool to annotate Mobile Element Insertions (MEIs). MEander takes VCFs from MEI callers, integrates disease-focused MEIs catalog, coding sequence coordinate, gnomAD-SV and generates reports for a quick overview of MEIs.

## Usage

### Installation
MEander requires some python dependencies, which can be installed through pip and conda.

``` bash
conda create -n MEander python=3.12
conda activate MEander

conda install -c bioconda bedtools -y
python -m pip install \
    polars==1.43.2 \
    fastexcel==0.20.2 \
    matplotlib==3.8.4 \
    dna-features-viewer==3.1.5 \
    pybedtools==0.12.0 \
    cyvcf2==0.32.1

```

Clone the repository
``` bash
git clone https://github.com/vincentxa847/MEander
python MEander.py --help
```

### Quick Start
``` bash
# NOTE: before running it, please see Input section to prepare the required inputs
INPUT="MELT.vcf_list.txt"
OUTPUT="sample"
OUTPUT_DIR="./content"

mkdir -p ${OUTPUT_DIR}

python MEander.py \
    --input ${INPUT} \
    --output ${OUTPUT} \
    --output-dir ${OUTPUT_DIR} \
    --mane ./MANE.gencode_v49.annotation.bed.gz \
    --pathogenic ./disease_focused_MEIs.xlsx \
    --clingen ./Clingen-Gene-Disease-Summary-2026-02-23.csv \
    --gnomad-constraint ./gnomad.v4.1.constraint_metrics.tsv \
    --control-dir ./gnomad/ \
    --gff ./gencode.v49.basic.annotation.gff3.gz \
    --mode CDS clingen MELT-QC \
    --window 20 \
    --position-to-map end \
    --caller MELT \
    --no-common
```

### All arguments
```
options:
  -h, --help            show this help message and exit
  -v, --version         Print version information

general options:
  --input INPUT         [REQUIRED] VCF or list of VCFs. If list, the file must end with .list.txt with one VCF path per line
  --output OUTPUT       [REQUIRED] Output prefix
  --output-dir OUTPUT_DIR
                        [REQUIRED] Directory for MEander outputs
  --verbose             Enable verbose/debug logging

mapping options:
  --position-to-map {start,end}
                        Position in the disease-associated MEI catalog (TSD start or end) to use. For MELT input, end is recommended (default: start)
  --mane MANE           [REQUIRED] MANE annotation BED file
  --pathogenic PATHOGENIC
                        [REQUIRED] Disease-focused MEI catalog file
  --clingen CLINGEN     [REQUIRED] ClinGen gene-disease summary file
  --gnomad-constraint GNOMAD_CONSTRAINT
                        [REQUIRED] gnomAD constraint metrics file
  --control-dir CONTROL_DIR
                        [REQUIRED] Directory with gnomAD resources
  --mode MODE [MODE ...]
                        Keep MEI in selected genomic features. [CDS, junction_donor, junction_acceptor, three_prime_UTR, five_prime_UTR]
  --clingen-only        Keep MEI in genes with ClinGen Gene-Disease Clinical Validity
  --qc {MELT,MEGAnE,xTea}
                        Remove low-confidence call. MELT: ASSESS >=3, split reads >2, FILTER=PASS; MEGAnE: FILTER=PASS. (default: none)
  --siteonly            When enabled, sample-level information is not collected
  --window WINDOW       Mapping window size (default: 20)
  --alt ALT             Field with MEI family (default: ALT)
  --clingen-classification CLASSIFICATION [CLASSIFICATION ...]
                        ClinGen classifications to include. [Definitive, Strong, Moderate, Limited, No Known Disease Relationship, Disputed, Refuted]

drawing options:
  --gff GFF             [REQUIRED] GFF3 annotation file
  --melt-meinfo         Handle MEI length and strand in MELT MEINFO explicitly
  --melen-field MELEN_FIELD
                        VCF field or INFO tag with MEI length
  --mestrand-field MESTRAND_FIELD
                        VCF field or INFO tag with MEI strand

report options:
  --no-common           Exclude common variants from the report
```

## Inputs
MEander takes 9 mandatory arguments as input:
1. INPUT : A VCF file or a list of VCF files (one VCF per line). MEander was tested on VCFs from MELT, MEGAnE and xTea.
2. OUTPUT : identifier for naming the output files.
3. OUTPUT_DIR : directory for the output files.
4. annotation resoueces
	- mapping
  		- [--pathogenic] Disease-focused MEI catalog file. Please download it from the manuscript.
      - [--clingen] Clingen Gene-Disease-Validity summary table from [ClinGen Gene-Disease Validity](https://search.clinicalgenome.org/kb/downloads#section_gene-disease-validity).
      - [--gnomad-constraint] gnomAD constraint table, download through `wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1.1/constraint/gnomad.v4.1.1.constraint_metrics.tsv.bgz`.
      - [--mane] BED file with MANE Select transcript exons with an additional 20 bp on each side of the exon–intron boundaries. See resource_meander.ipynb to prepare it.
      - [--control-dir] directory contains MEIs in gnomAD. See resource_meander.ipynb to prepare it.
	
    - drawing
    	- [--gff] GENCODE basic annotation, download through `wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gff3.gz`.

5. Options
	[--mode] : Mapping modes to run, which specify genomic feature(s) to retain when an MEI maps to them.

## Outputs
MEander produces these files in output directory:
- An HTML report for a quick overview of MEIs.
- A table with columns [source, gene,	feature, exon, mapped_pathogenic_MEI,	gnomad_AF, clingen annotation, lof.pLI] generated by the mapping module and all fields from the input VCF and INFO tags.
- A folder containing plots for each MEI, generated by the drawing module and named as `{CHROM}_{POS}_{ALT}.png`.
- A log file.