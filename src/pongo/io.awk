#' Read unique lines from a file.
#'
#' The lines are in arbitrary order.
#'
#' @param f File to read from. To ensure that lines are read from the beginning
#'   of the file, `f` will be closed before attempting to read from it.
#' @param arr Where to store the lines. `arr` will be cleared and then indexed
#'   from 1 to the number of lines read.
#' @param [skip] Regular expression used to skip lines.
#' @return Number of unique lines.
function read_lines_uniq(f, arr, skip,    old_rs, line, tmp, i) {
	old_rs = RS
	close(f)
	delete arr
	RS = "\n"

	while ((getline line < f) > 0) {
		if (line in tmp || (skip && line ~ skip)) {
			continue
		}

		arr[++i] = line
	}

	RS = old_rs

	return i
}
