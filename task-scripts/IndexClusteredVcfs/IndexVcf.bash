# Usage: IndexVcf.bash <VCF>
# Index a VCF using tabix

set -o errexit
set -o nounset
set -o pipefail

printf 'indexing %s\n' "$1" >&2
tabix "$1"
printf 'done\n' >&2
