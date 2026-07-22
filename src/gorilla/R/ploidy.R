# Read a sex ploidy file which is just the second and fifth columns of a PED
# file (sample ID, sex).
read_sex_ploidy <- function(path) {
    sex <- NULL
    sample_id <- NULL

    tmp <- fread(
        path,
        header = FALSE,
        sep = "\t",
        col.names = c("sample_id", "sex"),
        colClasses = c("character", "integer")
    )
    if (any(!tmp$sex %in% c(0L, 1L, 2L))) {
        stop("permitted values for sex are 0, 1, and 2")
    }

    tmp[sex == 0, sex := NA_integer_]
    tmp <- unique(tmp, by = "sample_id")
    setkey(tmp, sample_id)

    tmp
}
