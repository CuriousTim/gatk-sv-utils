# AddEnd2ToVcf.AddEnd2
# Args: <VCF> <OUTPUT_VCF>

set -o errexit
set -o nounset
set -o pipefail

in_vcf="$1"
out_vcf="$2"

# faster to skip bcftools parsing
bgzip -cd "${in_vcf}" \
  | awk -f /opt/task-scripts/utils/add_end2_to_vcf.awk - \
  | bgzip -c > "${out_vcf}"
bcftools index --tbi "${out_vcf}"
