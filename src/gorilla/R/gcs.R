get_gc_token <- function() {
    # GceToken does not implement validate, but a token on GCE can still expire
    # so we unconditionally refresh on GCE just to be safe
    if (is(gc_token, "GceToken") || !gc_token$validate()) {
        gc_token$refresh()
    }

    gc_token
}

gcs_download_file <- function(src, dest) {
    # bucket cannot contain forward slash
    m <- regexec("^gs://([^/]+)/(.+)$", src)
    if (m[[1]][[1]] == -1) {
        stop(sprintf("'%s' is not a valid Google Cloud Storage URI", src))
    }

    parts <- regmatches(src, m)[[1]]
    bucket <- parts[[2]]
    blob <- parts[[3]]

    req <- gargle::request_build(
        method = "GET",
        path = sprintf(
            "storage/v1/b/%s/o/%s",
            bucket,
            utils::URLencode(blob, reserved = TRUE)
        ),
        params = list(alt = "media"),
        token = get_gc_token(),
        base_url = "https://storage.googleapis.com"
    )

    res <- gargle::request_retry(
        req,
        httr::accept("application/octet-stream"),
        httr::write_disk(dest, overwrite = TRUE)
    )

    # don't use gargle::response_process
    # gargle enforces a JSON response, despite the underlying request
    # explicitly accepting something
    code <- httr::status_code(res)
    if (code < 200 || code >= 300) {
        stop(gargle::gargle_error_message(res))
    }

    invisible(dest)
}

get_gc_access_token <- function() {
    token <- get_gc_token()
    token$credentials$access_token
}
