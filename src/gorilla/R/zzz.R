.onLoad <- function(libname, pkgname) {
    env <- topenv()
    env$gc_token <- gargle::token_fetch(
        "https://www.googleapis.com/auth/cloud-platform"
    )
}
