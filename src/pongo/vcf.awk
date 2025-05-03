@namespace "vcf"

#' Parse a VCF meta-information or header line.
#'
#' @param line Line to parse.
#' @param a Array to fill with parsed information.
#' @return One of "unstructured", "structured", "header", or "unknown"
#'   indicating the type of the line.
function parse_header_line(line, a) {
	delete a
	if (line ~ /^##[^=]+=[^=]+$/) {
		_parse_unstructured(line, a)
		return "unstructured"
	}

	if (line ~ /^##[^=]+=<([^=<>]+=[^=<>]+)+>$/) {
		_parse_structured(line, a)
		return "structured"
	}

	if (line ~ /^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO/) {
		_parse_header(line, a)
		return "header"
	}

	return "unknown"
}

#' Parse VCF INFO field.
#'
#' @param info INFO field string.
#' @param a Array to fill with parsed information. The key-value pairs from the
#'   field will be accessible as key-values pairs in `a`.
function parse_info(info, a,    i, parts, kv) {
	delete a
	split(info, parts, /;/)
	for (i in parts) {
		n = split(parts[i], kv, /=/)
		if (n == 1) {
			a[kv[1]] = 1
		} else if (n == 2) {
			a[kv[1]] = kv[2]
		}
	}
}

#' Parse VCF genotype field.
#'
#' @param fmt FORMAT field for the record.
#' @param gt Genotype field for the sample.
#' @param a Array to populate with parsed information. The keys will be the
#'   ones from the FORMAT field and the values will be the corresponding values
#'   for the sample.
function parse_format(fmt, gt, a,    keys, values, i) {
	delete a
	split(fmt, keys, /:/)
	split(gt, values, /:/)
	for (i in keys) {
		a["fmt"][i] = keys[i]
		a["gt"][keys[i]] = values[i]
	}
}

#' Join VCF genotype parts.
#'
#' @param gt Array of genotype key-value pairs for the sample in the format
#'   created by `parse_format`.
#' @return Genotype string.
function make_format(gt,    fmt) {
	for (i in gt["fmt"]) {
		fmt = gt["gt"][gt["fmt"][i]] ":"
	}
	sub(/:$/, "", fmt)

	return fmt
}

function _parse_structured(x, a,    parts) {
	split(x, parts, /=/)
	sub(/^##/, "", parts[1])
	a[parts[1]] = parts[2]
}

function _parse_unstructured(x, a,    groups, parts, key, kv, pair) {
	match(x, /^##([^=]+)=<(.+)>$/, groups)
	key = groups[1]
	split(groups[2], parts, /,/)
	for (pair in parts) {
		split(pair, kv, /=/)
		a[key][kv[1]] = kv[2]
	}
}

function _parse_header(x, a,    parts, n, i) {
	n = split(x, parts, /\t/)
	for (i = 10; i <= n; ++i) {
		a["samples"]["pos"][i] = parts[i]
		a["samples"]["id"][parts[i]] = i
	}
}
