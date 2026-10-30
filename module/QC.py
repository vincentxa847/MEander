'''
QC module

- MELT (criteria from https://academic.oup.com/nar/article/50/5/2493/6536893)
    - ASSESS score ≥3 (at least one-side TSD evidence)
    - split reads >2 (at least two split-read support)
    - FILTERS=PASS
- MEGAnE (criteria from https://www.nature.com/articles/s41588-023-01390-2)
    - FILTER=PASS
- xTea (criteria from https://link.springer.com/article/10.1186/s13100-021-00256-w)
    - clipped reads ≥3
    - discordant pairs ≥5

INPUT  : cyvcf2 variant record
OUTPUT : Status of PASS or NOT PASS QC 
'''


def passes_qc(v, qc):
    """
    Return True if a Variant passes QC.
    """

    if "MELT" in qc:

        # cyvcf2: both PASS and '.' are represented as None
        if v.FILTER is not None:
            return False

        sr = v.INFO.get("SR", 0)
        assess = v.INFO.get("ASSESS", 0)

        if sr <= 2 or assess < 3:
            return False

    if "MEGAnE" in qc:

        # cyvcf2: both PASS and '.' are represented as None
        if v.FILTER is not None:
            return False

    if "xTea" in qc:

        lclip = v.INFO.get("LCLIP", 0)
        rclip = v.INFO.get("RCLIP", 0)

        ldisc = v.INFO.get("LDISC", 0)
        rdisc = v.INFO.get("RDISC", 0)

        clip = lclip + rclip
        disc = ldisc + rdisc

        if clip < 3 or disc < 5:
            return False

    return True
