parse_pkgs <- function(x) {
    strsplit(x, split = "@", fixed = TRUE)
}

nproc <- function() {
    withCallingHandlers(as.integer(system2("nproc", stdout = TRUE)),
                        warning = function(w) 1L)
}

NCPUS <- nproc()

install <- function(x) {
    withCallingHandlers(
        remotes::install_version(
            x[[1]],
            version = x[[2]],
            upgrade = "always",
            repos = c("https://cloud.r-project.org", "https://bioconductor.org/packages/3.22/bioc"),
            INSTALL_opts = c("--no-docs", "--no-html", "--no-data", "--no-help",
                             "--no-demo", "--without-keep.source"),
            Ncpus = NCPUS
        ),
        warning = function(w) stop(w)
    )
}

# Packages --------------------------------------------------------------------
cran_pkgs <- c("data.table@1.17.8",
               "R.utils@2.13.0",
               "matrixStats@1.5.0") |>
  parse_pkgs()
bioc_pkgs <- c("Rsamtools@2.26.0",
               "GenomicRanges@1.62.0",
	       "rtracklayer@1.70.0") |>
  parse_pkgs()

invisible(lapply(cran_pkgs, install))
invisible(lapply(bioc_pkgs, install))
