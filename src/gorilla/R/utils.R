TABIX_MAX_SEQLEN <- 536870912L

#' Expand ranges by a fraction of their size.
#'
#' The ranges are expanded to a minimum start of 1 and a maximum end of
#' 536,870,912 which are the limits of the tabix library.
#'
#' @param x `data.table` The ranges to expand. It must have at least the
#'   columns `chr`, `start`, and `end`.
#' @param f `double` The fraction by which to expand the ranges.
#' @returns `data.table` The expanded ranges.
#' @export
expand_tabix_ranges <- function(x, f) {
    start <- NULL
    end <- NULL

    sizes <- x[, end - start + 1L]
    pad <- ceiling(sizes * f)
    x[, start := pmax(1L, start - pad)]
    x[, end := pmin(end + pad, TABIX_MAX_SEQLEN)]

    x
}
