# Usage: IndexVcf.bash <VCF>
# Index a VCF using tabix

set -o errexit
set -o nounset
set -o pipefail

vcf="$1"
vcf_bn="$(basename "${vcf}")"

# The index file will be placed in the current directory so that it will be
# delocalized to the top level task output directory.
printf 'indexing %s\n' "${vcf}" >&2
tabix "${vcf}"
mv "${vcf}.tbi" "${vcf_bn}.tbi"
printf 'done\n' >&2
