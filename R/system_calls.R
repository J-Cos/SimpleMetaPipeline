#' Run an external command with error checking
#'
#' A wrapper around [system2()] that checks the exit code and raises an
#' informative error on failure. All external tool calls in the pipeline go
#' through this function.
#'
#' @param command The command to run (e.g. `"vsearch"`).
#' @param args Character vector of arguments.
#' @param step_name A human-readable label for log messages.
#' @param stdout Passed to [system2()]. Default `""` (print to console).
#' @param stderr Passed to [system2()]. Default `""` (print to console).
#'
#' @return Invisible character vector of stdout lines (if captured).
#' @keywords internal
smp_system <- function(command, args, step_name, stdout = "", stderr = "") {
    cli::cli_alert("Running {.val {command}} for {.val {step_name}}")

    result <- system2(command, args, stdout = stdout, stderr = stderr)

    # system2 returns the exit code when stdout/stderr are ""
    # and returns output with an attribute "status" when capturing
    if (is.character(result)) {
        exit_code <- attr(result, "status")
    } else {
        exit_code <- result
        result <- character(0)
    }

    if (!is.null(exit_code) && exit_code != 0) {
        cli::cli_abort(c(
            "{.val {command}} failed with exit code {exit_code} at step {.val {step_name}}",
            "i" = if (length(result) > 0) {
                paste("Last output:", paste(utils::tail(result, 5), collapse = "\n"))
            } else {
                "No captured output. Check console for stderr."
            }
        ))
    }

    invisible(result)
}


#' Run an external command via system() with error checking
#'
#' For commands that need the full shell (e.g. piping, `&&`). Prefer
#' [smp_system()] when possible.
#'
#' @param command_string The full shell command string.
#' @param step_name A human-readable label for log messages.
#'
#' @return Invisible exit code (0).
#' @keywords internal
smp_shell <- function(command_string, step_name) {
    cli::cli_alert("Running shell command for {.val {step_name}}")
    exit_code <- system(command_string)

    if (exit_code != 0) {
        cli::cli_abort(c(
            "Shell command failed with exit code {exit_code} at step {.val {step_name}}",
            "i" = "Command: {command_string}"
        ))
    }

    invisible(exit_code)
}
