@namespace "logging"

#' Return the current time in UTC as a formatted string.
function datetime() {
	return awk::strftime("%FT%H:%M:%S%Z", awk::systime(), 1)
}

#' Log an INFO level message.
#'
#' @param msg Message to log.
#' @param [file] Where to write the message. Default is stderr.
function log_info(msg, file) {
	_log_msg(msg, "INFO", file)
}

#' Log a WARN level message.
#'
#' @param msg Message to log.
#' @param [file] Where to write the message. Default is stderr.
function log_warn(msg, file) {
	_log_msg(msg, "WARN", file)
}

#' Log an ERROR level message.
#'
#' @param msg Message to log.
#' @param [file] Where to write the message. Default is stderr.
function log_err(msg, file) {
	_log_msg(msg, "ERROR", file)
}

function _log_msg(msg, level, file) {
	printf "[%s %s] %s\n", datetime(), level, msg > (file ? file : "/dev/stderr")
}
