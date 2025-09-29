"""
Add FILTERs to a VCF

usage: python add_vcf_filters.py <invcf> <outvcf> <filters> [<new_headers>]

<invcf>          VCF to update.
<outvcf>         Where to write the updated VCF.
<filters>        A TSV with variant ID in the first column and semicolon-separated
                 FILTERs in the second column. No header.
[<new_headers>]  Additional header lines to add. Required if FILTERs not already
                 defined in <invcf> are being added.

There is no error checking so get the inputs right.
"""

import sys

import pysam
from pysam import VariantFile


def load_filters(path):
    filters = dict()
    with open(path, mode="r", encoding="utf-8") as f:
        for line in f:
            vid, fil = line.rstrip().split("\t")
            old_fil = filters.get(vid, set())
            filters[vid] = old_fil.union(set(fil.split(";")))

    return filters


def main():
    invcf = VariantFile(sys.argv[1], mode="r")
    filters = load_filters(sys.argv[3])
    header = invcf.header
    if len(sys.argv) == 5:
        with open(sys.argv[4], mode = "r", encoding = "utf-8") as f:
            for line in f:
                header.add_line(line.rstrip())
    outvcf = VariantFile(sys.argv[2], mode="w", header=header)

    try:
        for rec in invcf.fetch():
            if rec.id in filters:
                fil = filters[rec.id]
                for x in fil:
                    rec.filter.add(x)
            outvcf.write(rec)
    finally:
        invcf.close()
        outvcf.close()

    pysam.tabix_index(sys.argv[2], force = True, preset = "vcf")


if __name__ == "__main__":
    main()
