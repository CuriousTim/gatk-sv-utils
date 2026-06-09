# Some functions to help use the system tabix to query tabix-indexed file. It
# is easier to use tabix directly for accessing files hosted on cloud storage
# than it is to use Rsamtools, which uses the htslib library bundled with
# Rhtslib. It is important that the system tabix uses a version of htslib with
# the fixes for silently dropping records from remote files, which were
# completed in commit e805d4031a17b38798d64cae4aaeae1f3cff2be6.

#' @export
query.tabix_handle <- function(x, contig, start, end) {
    # tabix will save a remote index file to the current directory so we need
    # to change the current directory to the cache directory
    old_wd <- setwd(x$cachedir)
    on.exit(setwd(old_wd), add = TRUE)

    results <- tempfile(tmpdir = ".")
    rc <- system2(
        "tabix",
        args = c(
            "--print-header",
            x$path,
            sprintf("%s:%d-%d", contig, start, end)
        ),
        stdout = results,
        env = sprintf("GCS_OAUTH_TOKEN=%s", get_gc_access_token())
    )

    if (rc != 0) {
        stop("tabix failed")
    }

    file.path(x$cachedir, results)
}

new_tabix_handle <- function(path, cachedir) {
    if (is.null(cachedir)) {
        cachedir <- tempfile("tabix", tempdir(TRUE))
    }

    tryCatch(
        mkdir(cachedir),
        error = function(e) {
            stop("failed to create tabix cache directory")
        }
    )

    structure(list(path = path, cachedir = cachedir), class = "tabix_handle")
}
