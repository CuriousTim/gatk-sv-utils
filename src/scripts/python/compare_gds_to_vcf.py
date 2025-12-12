"""
Compare manually reviewed genomic disoder CNVs against a VCF.

usage: compare_gds_to_vcf.py <tp_gds> <vcf> <matches> <mismatches>

<gds>         Table of genomic disorders with carriers and non-carriers.
<vcf>         VCF.
<matches>     Variants in VCF for which all carriers match the
              carriers of the GD.
<mismatches>  Variant in VCF for which one or more carriers do
              not match the carriers of the GD.

The format of <gds> is a tab-delimited file with columns:
chr
start (1-based)
end (closed end)
GD ID
cluster ID
SV type
Manually validated carrriers (comma-separated)
Manually validated non-carriers (comma-separated)

<matches> will contain the IDs of variants in the VCF that match a GD CNV and
for which the set of the GD carriers is equal to the set of VCF variant
carriers.

<mismatches> will contain a table of VCF variants that match a GD CNV, but
there is discordance between the VCF variant carriers and the GD CNV
carriers.
The format of <mismatches> is a tab-delimited file with columns:
VCF ID
VCF start
VCF end
SV type
GD ID
VCF carriers that are manually reviewed GD carriers (comma-separated)
VCF carriers that are manually reviewed non-GD carriers (comma-separated)
VCF carriers that are not in the previous two categories (comma-separated)
"""

import sys
import enum
from enum import Enum

import pysam
from pysam import VariantFile


class SVType(Enum):
    DEL = enum.auto()
    DUP = enum.auto()
    INS = enum.auto()
    CTX = enum.auto()
    CPX = enum.auto()
    INV = enum.auto()
    BND = enum.auto()

    @classmethod
    def from_string(cls, s):
        if not isinstance(s, str):
            raise TypeError("`s` must be a string")
        e = cls.__members__.get(s)
        if e is None:
            raise ValueError(f"{s} is not a valid SV type string")
        return e

    def __str__(self):
        return f"{self.name}"


class Interval:
    def __init__(self, start, end):
        Interval._validate(start, end)
        self._start = start
        self._end = end

    @property
    def start(self):
        return self._start

    @start.setter
    def start(self, value):
        if not Interval._is_positive_int(value):
            raise ValueError("start must be a positive integer")
        Interval._assert_well_ordered(value, self._end)

        self._start = value

    @property
    def end(self):
        return self._end

    @end.setter
    def end(self, value):
        if not Interval._is_positive_int(value):
            raise ValueError("end must be a positive integer")
        Interval._assert_well_ordered(self._start, value)

        self._end = value

    @staticmethod
    def _validate(start, end):
        if not Interval._is_positive_int(start) or not Interval._is_positive_int(end):
            raise ValueError("start and end must be a positive integers")
        Interval._assert_well_ordered(start, end)

    @staticmethod
    def _assert_well_ordered(start, end):
        if start > end:
            raise ValueError("start must be less than or equal to end")

    def __str__(self):
        return f"{self.start}\t{self.end}"

    def overlap(self, other):
        return max(0, min(self._end, other._end) - max(self._start, other._start) + 1)

    def size(self):
        return self._end - self._start + 1

    @staticmethod
    def _is_positive_int(x):
        return x and isinstance(x, int) and x > 0


class CNV:
    def __init__(self, contig, interval, svtype):
        if not CNV._is_valid_svtype(svtype):
            raise ValueError("Only DEL and DUP are permitted")
        self.contig = contig
        self.interval = interval
        self.svtype = svtype

    @classmethod
    def from_parts(cls, contig, start, end, svtype):
        return cls(contig, Interval(start, end), SVType.from_string(svtype))

    def matches(self, other, min_ovp=0.5):
        if self.contig != other.contig or self.svtype is not other.svtype:
            return False

        ovp = self.interval.overlap(other.interval)

        return (
            ovp / self.interval.size() >= min_ovp
            and ovp / other.interval.size() >= min_ovp
        )

    def __str__(self):
        return f"{self.contig}\t{self.interval}\t{self.svtype}"

    @staticmethod
    def _is_valid_svtype(x):
        return isinstance(x, SVType) and (x is SVType.DEL or x is SVType.DUP)


class VCFRecord:
    def __init__(self, vr):
        self.variant = CNV.from_parts(vr.contig, vr.pos, vr.stop, vr.info["SVTYPE"])
        self.vid = vr.id
        self.carriers = set(VCFRecord._gt_to_carriers(vr.samples))

    def overlaps_gd(self, gdrecord, min_ovp):
        return self.variant.matches(gdrecord.variant, min_ovp)

    @staticmethod
    def _gt_to_carriers(x):
        return [
            sample for sample in x.keys() if VCFRecord._gt_has_allele(x[sample]["GT"])
        ]

    @staticmethod
    def _gt_has_allele(gt):
        return (
            len(gt) == 2
            and (gt[0] is not None)
            and (gt[1] is not None)
            and (gt[0] == 1 or gt[1] == 1)
        )


class GDRecord:
    def __init__(
        self, contig, start, end, gdid, clustid, svtype, carriers, non_carriers
    ):
        self.variant = CNV.from_parts(contig, int(start), int(end), svtype)
        self.gdid = gdid
        self.clustid = None
        self.carriers = set(carriers.split(","))
        self.non_carriers = set(non_carriers.split(","))


class GDTable:
    def __init__(self, records, clusters):
        self.records = records
        self.clusters = clusters

    @classmethod
    def from_file(cls, path):
        records = list()
        clusters = dict()
        with open(path, mode="r", encoding="utf-8") as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                rec = GDRecord(*fields)
                records.append(rec)
                if rec.clustid is not None:
                    clusters.get(rec.clustid, set()).update(rec.carriers)
        return cls(records, clusters)

    def iter_records(self):
        for rec in self.records:
            yield rec


class GDComparator:
    def __init__(self, gdtable, variantfile):
        self.gdtable = gdtable
        self.variantfile = variantfile
        self.variantfile_samples = set(variantfile.header.samples)

    def get_matches(self, min_ovp):
        for gdrec in self.gdtable.iter_records():
            gdvar = gdrec.variant
            if not gdvar.contig in self.variantfile.index:
                continue
            for vfilerec in self.variantfile.fetch(
                gdvar.contig, gdvar.interval.start, gdvar.interval.end
            ):
                if not (
                    "PASS" in vfilerec.filter or "FAIL_MANUAL_REVIEW" in vfilerec.filter
                ):
                    continue
                try:
                    vcfrec = VCFRecord(vfilerec)
                except ValueError:
                    # raised by CNV.__init__() so assume SV type is not CNV, which might not be the problem
                    continue
                if not vcfrec.overlaps_gd(gdrec, min_ovp):
                    continue

                c1, c2, c3, c4 = GDComparator._get_sample_overlaps(vcfrec, gdrec, self.gdtable)
                yield (vcfrec, gdrec, c1, c2, c3, c4)

    @staticmethod
    def _get_sample_overlaps(vcfrec, gdrec, gdtable):
        if gdrec.clustid is None:
            cluster_carriers = set()
        else:
            # every cluster ID should have an entry
            cluster_carriers = gdtable.clusters[gdrec.clustid]
        # VCF carriers that are true GD carriers
        c1 = vcfrec.carriers.intersection(gdrec.carriers)
        # VCF carriers that are true GD non-carriers for this GD, but a carrier
        # for a different GD in the same cluster
        c2 = vcfrec.carriers.intersection(gdrec.non_carriers).intersection(
            cluster_carriers
        )
        # VCF carriers that are true GD non-carriers for this GD and not a
        # carrier for a different GD in the same cluster
        c3 = vcfrec.carriers.intersection(gdrec.non_carriers).difference(
            cluster_carriers
        )
        # VCF carriers that are not otherwise classified
        c4 = vcfrec.carriers.difference(
            gdrec.carriers.intersection(gdrec.non_carriers)
        )
        return (c1, c2, c3, c4)


def compare(gds, vcf, matches, mismatches):
    gdc = GDComparator(gds, vcf)
    rm_vids = set()
    other_vids = set()
    with (
        open(matches, mode="w", encoding="utf-8") as f,
        open(mismatches, mode="w", encoding="utf-8") as g,
    ):
        for vcfrec, gdrec, c1, c2, c3, c4 in gdc.get_matches(0.5):
            if len(c3) == 0 and len(c4) == 0:
                rm_vids.add(vcfrec.vid)
            else:
                other_vids.add(vcfrec.vid)
                g.write(
                    f"{vcfrec.vid}\t{vcfrec.variant}\t{gdrec.gdid}\t{','.join(c1)}\t{','.join(c2)}\t{','.join(c3)}\t{','.join(c4)}\n"
                )
        rm_vids.difference_update(other_vids)
        f.writelines(rm_vids)


def main():
    gds = GDTable.from_file(sys.argv[1])
    vcf = VariantFile(sys.argv[2], mode="r")
    matches = sys.argv[3]
    mismatches = sys.argv[4]

    compare(gds, vcf, matches, mismatches)


if __name__ == "__main__":
    main()
