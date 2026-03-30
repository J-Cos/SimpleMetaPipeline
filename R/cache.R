#' Cache a pipeline step result to disk
#'
#' Saves the result of a pipeline step as an RDS file in the
#' IntermediateOutputs directory, enabling resumption of interrupted runs.
#'
#' @param result The result object to cache.
#' @param config An `smp_config` object.
#' @param step_name A short identifier for this step (e.g. `"dada2"`).
#' @return Invisible `result`.
#' @keywords internal
smp_cache <- function(result, config, step_name) {
    out_path <- file.path(
        config$path, "IntermediateOutputs",
        paste0(config$dataname, "_", step_name, ".RDS")
    )
    saveRDS(result, out_path)
    cli::cli_alert_info("Cached {.val {step_name}} to {.path {out_path}}")
    invisible(result)
}


#' Load all cached intermediate results
#'
#' Scans the IntermediateOutputs directory for previously cached results and
#' returns them in a named list. Used by [smp_run()] to resume from the last
#' completed step. Supports both new (`_dada2.RDS`) and legacy
#' (`_DadaOutput.RDS`) naming conventions for backward compatibility.
#'
#' @param config An `smp_config` object.
#' @return A named list with one entry per cached step. Missing steps are
#'   `NULL`.
#' @keywords internal
smp_load_cache <- function(config) {
    # New names → legacy (old pipeline) names
    step_aliases <- list(
        cutadapt = c("cutadapt", "CutadaptOutput"),
        dada2    = c("dada2", "DadaOutput"),
        lulu1    = c("lulu1", "LuluOutput1"),
        cluster  = c("cluster", "ClusterOutput"),
        lulu2    = c("lulu2", "LuluOutput2"),
        idtaxa   = c("idtaxa", "IdtaxaOutput"),
        blast    = c("blast", "BlastOutput")
    )

    cache <- setNames(vector("list", length(step_aliases)), names(step_aliases))

    for (step in names(step_aliases)) {
        for (alias in step_aliases[[step]]) {
            rds_path <- file.path(
                config$path, "IntermediateOutputs",
                paste0(config$dataname, "_", alias, ".RDS")
            )
            if (file.exists(rds_path)) {
                cache[[step]] <- readRDS(rds_path)
                if (alias != step) {
                    cli::cli_alert_info(
                        "Loaded legacy cached {.val {alias}} as {.val {step}}"
                    )
                } else {
                    cli::cli_alert_info(
                        "Loaded cached {.val {step}} from previous run"
                    )
                }
                break
            }
        }
    }

    cache
}


#' Clear cached intermediate results
#'
#' Removes all cached RDS files for a given dataset from the
#' IntermediateOutputs directory. Use this before re-running the pipeline on
#' the same dataset.
#'
#' @param config An `smp_config` object.
#' @return Invisible `TRUE`.
#' @export
smp_clear_cache <- function(config) {
    pattern <- paste0("^", config$dataname, "_.*\\.RDS$")
    cache_dir <- file.path(config$path, "IntermediateOutputs")
    files <- list.files(cache_dir, pattern = pattern, full.names = TRUE)
    if (length(files) > 0) {
        file.remove(files)
        cli::cli_alert_success(
            "Removed {length(files)} cached file{?s} for {.val {config$dataname}}"
        )
    } else {
        cli::cli_alert_info("No cached files found for {.val {config$dataname}}")
    }
    invisible(TRUE)
}
