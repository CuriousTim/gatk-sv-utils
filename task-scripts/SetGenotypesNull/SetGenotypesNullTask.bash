# SetGenotypesNull.SetGenotypesNullTask
# Args: <VCF> <OUTPUT_VCF> <SAMPLES_LIST> <CONTIGS>

set -o errexit
set -o nounset
set -o pipefail

in_vcf="$1"
out_vcf="$2"
samples="$3"
contigs="$4"

if [[ $(wc -l "${samples}" | awk '{print $1}') -eq 0 ]]; then
  cp "${in_vcf}" "${out_vcf}"
else
  bgzip -cd "${in_vcf}" \
    gawk -f /opt/task-scripts/utils/set_gt_null.awk "${contigs}" "${samples}" - \
    bgzip -c > "${out_vcf}"
fi
bcftools index --tbi "${out_vcf}"
