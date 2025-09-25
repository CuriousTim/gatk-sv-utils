"""
Nullify sample genotypes driven by batch effects.

usage: nullify_batch_effect_gts.py <invcf> <outvcf> <batches> <variants>

<invcf>     VCF to update.
<outvcf>    Where to write the updated VCF.
<batches>   A TSV with sample ID in first column and batch ID in second column.
            No header.
<variants>  A TSV with variant ID in first column and semicolon separated batch
            IDs in second column. No header.

There is no error checking so get the inputs right.
"""

import sys

import pysam
from pysam import VariantFile


def load_batch_table(path, samples_to_keep):
    batches = dict()
    with open(path, mode="r", encoding="utf-8") as f:
        for line in f:
            sid, bid = line.rstrip().split("\t")
            if sid in samples_to_keep:
                samples = batches.get(bid, set())
                samples.add(sid)
                batches[bid] = samples

    return batches


def load_variant_table(path):
    vids = dict()
    with open(path, mode="r", encoding="utf-8") as f:
        for line in f:
            vid, bids = line.rstrip().split("\t")
            vids[vid] = bids.split(";")

    return vids


def main():
    invcf = VariantFile(sys.argv[1], mode="r")
    outvcf = VariantFile(sys.argv[2], mode="w", header=invcf.header)
    batch_map = load_batch_table(sys.argv[3], invcf.header.samples)
    vid_map = load_variant_table(sys.argv[4])

    try:
        for rec in invcf.fetch():
            if rec.id in vid_map:
                bids = vid_map[rec.id]
                for bid in bids:
                    for sid in batch_map[bid]:
                        if rec.info["SVTYPE"] == "CNV":
                            rec.samples[sid]["CN"] = None
                        else:
                            rec.samples[sid]["GT"] = (None, None)
            outvcf.write(rec)
    finally:
        invcf.close()
        outvcf.close()

    pysam.tabix_index(sys.argv[2], preset = "vcf")


if __name__ == "__main__":
    main()
