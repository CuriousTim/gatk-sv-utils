# SplitVcfByContig.GetContigFromVcf
# Args: <VCF> <CONTIG> <OUTPUT_VCF>

set -o errexit
set -o nounset
set -o pipefail

in_vcf="$1"
contig="$2"
out_vcf="$3"

bcftools index --stats "${in_vcf}" \
  | awk -F'\t' '$1 == a{b=1} END{if(!b){printf "%s not in VCF\n", a > "/dev/stderr"; exit 1}}' a="${contig}"
bcftools head "${in_vcf}" \
  | awk '
    /^##contig=<ID=/ {
      match($0, /<ID=[^,>]+/)
      contig = substr($0, RSTART + 4, RLENGTH - 4)
      if (contig == target) {
        print
      }
      next
    }
    1' target="${contig}" > 'new_header.txt'
bcftools view --output-type u --regions "${contig}" "${in_vcf}" \
  | bcftools reheader --header 'new_header.txt' - \
  | bcftools view --output-type z --output "${out_vcf}" -
bcftools index --tbi "${out_vcf}"
