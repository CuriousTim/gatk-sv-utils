"""
Add INFOs to a VCF

usage: python add_vcf_infos.py <invcf> <outvcf> <infos> [<new_headers>]

<invcf>          VCF to update.
<outvcf>         Where to write the updated VCF.
<infos>          A TSV with variant ID in the first column, the INFO key in the
                 second column, and the INFO value in the thirds. If a flag is
                 being added, the third column of the row can be omitted.
[<new_headers>]  Additional header lines to add. Required if INFOs not already
                 defined in <invcf> are being added.

There is no error checking so get the inputs right.
"""

import sys

import pysam
from pysam import VariantFile


def load_infos(path):
    infos = dict()
    with open(path, mode="r", encoding="utf-8") as f:
        for line in f:
            fields = line.rstrip().split("\t")
            vid = infos.get(fields[0], dict())
            vid[fields[1]] = True if len(fields) == 2 else fields[2]
            infos[fields[0]] = vid

    return infos


def main():
    invcf = VariantFile(sys.argv[1], mode="r")
    infos = load_infos(sys.argv[3])
    header = invcf.header
    if len(sys.argv) == 5:
        with open(sys.argv[4], mode = "r", encoding = "utf-8") as f:
            for line in f:
                header.add_line(line.rstrip())
    outvcf = VariantFile(sys.argv[2], mode="w", header=header)

    try:
        for rec in invcf.fetch():
            rec.info.update(infos.get(rec.id, dict()))
            outvcf.write(rec)
    finally:
        invcf.close()
        outvcf.close()

    pysam.tabix_index(sys.argv[2], force = True, preset = "vcf")


if __name__ == "__main__":
    main()
