# SplitVcfByContig.GetContigFromVcf
# Args: <VCF> <CONTIG> <OUTPUT_VCF>

set -o errexit
set -o nounset
set -o pipefail

bcftools view --output "$3" --output-type z --regions "$2" "$1"
bcftools index --tbi "$3"
