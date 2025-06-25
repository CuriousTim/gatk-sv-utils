parse_pkgs <- function(x) {
    strsplit(x, split = "@", fixed = TRUE)
}

nproc <- function() {
    withCallingHandlers(as.integer(system2("nproc", stdout = TRUE)),
                        warning = function(w) 1L)
}

NCPUS <- nproc()

install <- function(x) {
    remotes::install_version(
        x[[1]],
        version = x[[2]],
        upgrade = "always",
        repos = c("https://cloud.r-project.org", "https://bioconductor.org/packages/3.21"),
        INSTALL_opts = c("--no-docs", "--no-html", "--no-data", "--no-help",
                         "--no-demo", "--without-keep.source"),
        Ncpus = NCPUS
    )
}

# Packages --------------------------------------------------------------------
cran_pkgs <- c("data.table@1.17.6",
               "optparse@1.7.5") |>
  parse_pkgs()
bioc_pkgs <- c("Rsamtools@2.24.0",
               "GenomicRanges@1.60.0") |>
  parse_pkgs()

invisible(lapply(cran_pkgs, install))
invisible(lapply(bioc_pkgs, install))
