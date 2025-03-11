# MoveVcfFiltersToInfo.MoveVcfFilters
# Args: <VCF> <OUTPUT_VCF>

set -o errexit
set -o nounset
set -o pipefail

in_vcf="$1"
out_vcf="$2"

bgzip -cd "${in_vcf}" \
  | gawk -f /opt/task-scripts/utils/move_vcf_filters.awk - \
  | bgzip -c > "${out_vcf}"
bcftools index --tbi "${out_vcf}"
