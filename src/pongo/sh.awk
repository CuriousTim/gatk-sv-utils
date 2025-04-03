@namespace "sh"

#' Quote a string for use with the shell.
#'
#' All single quotes in the string will be escaped and then the entire string
#' will be wrapped in single quotes. Existing single quote characters in the
#' string are escaped by replacing each one with the sequence '\'' (all four
#' characters are literal). The first quote closes the previous quoted string.
#' The backslash plus quote escapes the quote and the shell will remove the
#' backslash. The final quote starts a new quoted string. After replacement and
#' wrapping everything in quotes, the resulting string will be a sequence of
#' single-quoted strings with intervening escaped single quotes, which the
#' shell will concatenate into a single string. All characters in single-quotes
#' are interpreted literally by the shell.
#'
#' @param x String to quote.
function quote(x) {
	gsub(/'/, "'\\''", x)
	return "'" x "'"
}
