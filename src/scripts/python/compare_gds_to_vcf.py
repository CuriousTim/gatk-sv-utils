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


class GDFileRecord:
    def __init__(self, contig, start, end, gdid, svtype, carriers, non_carriers):
        self.variant = CNV.from_parts(contig, int(start), int(end), svtype)
        self.gdid = gdid
        self.carriers = set(carriers.split(","))
        self.non_carriers = set(non_carriers.split(","))


class GDFile:
    def __init__(self, path):
        self.path = path
        self.fp = None

    def __enter__(self):
        self.fp = open(self.path, encoding="utf-8")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.fp:
            self.fp.close()
            self.fp = None

    def records(self):
        if self.fp and not self.fp.closed:
            for line in self.fp:
                fields = line.rstrip("\n").split("\t")
                yield GDFileRecord(*fields)


class GDComparator:
    def __init__(self, gdfile, variantfile):
        self.gdfile = gdfile
        self.variantfile = variantfile
        self.variantfile_samples = set(variantfile.header.samples)

    def get_matches(self, min_ovp):
        with self.gdfile as g:
            for gdrecord in g.records():
                gdvar = gdrecord.variant
                if not gdvar.contig in self.variantfile.index:
                    continue
                # just in case there are samples in the GD file that aren't in the VCF
                gdrecord.carriers = gdrecord.carriers.intersection(
                    self.variantfile_samples
                )
                gdrecord.carriers = gdrecord.non_carriers.intersection(
                    self.variantfile_samples
                )
                for rec in self.variantfile.fetch(
                    gdvar.contig, gdvar.interval.start, gdvar.interval.end
                ):
                    if not ("PASS" in rec.filter or "FAIL_MANUAL_REVIEW" in rec.filter):
                        continue
                    try:
                        vcfrecord = VCFRecord(rec)
                    except ValueError:
                        # raised by CNV.__init__() so assume SV type is not CNV, which might not be the problem
                        continue
                    if vcfrecord.overlaps_gd(gdrecord, min_ovp):
                        yield (gdrecord, vcfrecord)


def compare(gds, vcf, matches, mismatches):
    """Compare genomic disoder CNVs to variants in a VCF.

    For each CNV in `gds`, find each matching CNV in `vcf` and compare the
    carriers of the two. If the two have the same set of carriers, write the VCF
    ID to `matches`. If the two do not have the same set of carriers, write the
    variant and the discordant carriers to `mismatches`.
    """
    gdc = GDComparator(gds, vcf)
    with (
        open(matches, mode="w", encoding="utf-8") as f,
        open(mismatches, mode="w", encoding="utf-8") as g,
    ):
        for gdrec, vcfrec in gdc.get_matches(0.5):
            # GD CNV carriers that are VCF CNV carriers
            common_tp = gdrec.carriers.intersection(vcfrec.carriers)
            # GD CNV non-carriers that are VCF CNV carriers
            common_tn = gdrec.non_carriers.intersection(vcfrec.carriers)
            # VCF CNV carriers not in previous sets
            others = vcfrec.carriers.difference(common_tp.union(common_tn))
            if len(common_tn) == 0 and len(others) == 0:
                f.write(f"{vcfrec.vid}\n")
            else:
                g.write(
                    f"{vcfrec.vid}\t{vcfrec.variant}\t{gdrec.gdid}\t{",".join(common_tp)}\t{",".join(common_tn)}\t{",".join(others)}\n"
                )


def main():
    gds = GDFile(sys.argv[1])
    vcf = VariantFile(sys.argv[2], mode="r")
    matches = sys.argv[3]
    mismatches = sys.argv[4]

    compare(gds, vcf, matches, mismatches)


if __name__ == "__main__":
    main()
