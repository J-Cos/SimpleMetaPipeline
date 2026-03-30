#' Validate a pipeline configuration before running
#'
#' Checks that all required files, directories, and external tools exist
#' before committing to a potentially long pipeline run. Called automatically
#' by [smp_run()].
#'
#' @param config An `smp_config` object.
#' @return Invisible `TRUE` if all checks pass.
#' @export
smp_validate <- function(config) {
    errors <- character()
    warnings <- character()

    # --- Directory structure ---
    fastq_dir <- file.path(config$data_path, "FASTQs", config$dataname)
    if (!dir.exists(fastq_dir)) {
        errors <- c(errors, paste0(
            "FASTQ directory not found: ", fastq_dir
        ))
    } else {
        run_dirs <- list.files(fastq_dir)
        if (length(run_dirs) < config$n_runs) {
            errors <- c(errors, paste0(
                "Expected ", config$n_runs, " run directories in ", fastq_dir,
                " but found ", length(run_dirs), ": ",
                paste(run_dirs, collapse = ", ")
            ))
        }
    }

    int_out_dir <- file.path(config$path, "IntermediateOutputs")
    if (!dir.exists(int_out_dir)) {
        errors <- c(errors, paste0(
            "IntermediateOutputs directory not found: ", int_out_dir
        ))
    }

    results_dir <- file.path(config$path, "Results")
    if (!dir.exists(results_dir)) {
        errors <- c(errors, paste0(
            "Results directory not found: ", results_dir
        ))
    }

    # --- Parameter consistency ---
    if (length(config$trunc_len) != config$n_runs) {
        errors <- c(errors, paste0(
            "trunc_len has ", length(config$trunc_len),
            " entries but n_runs = ", config$n_runs
        ))
    }

    if (length(config$trim_left) != config$n_runs) {
        errors <- c(errors, paste0(
            "trim_left has ", length(config$trim_left),
            " entries but n_runs = ", config$n_runs
        ))
    }

    if (config$read_type == "single-read") {
        for (i in seq_along(config$trunc_len)) {
            if (length(config$trunc_len[[i]]) > 1) {
                warnings <- c(warnings, paste0(
                    "read_type is 'single-read' but trunc_len[[", i,
                    "]] has ", length(config$trunc_len[[i]]),
                    " values. Only the first will be used."
                ))
            }
        }
    }

    # --- Classifier / BLAST DB ---
    if (config$use_idtaxa) {
        if (is.null(config$training_set)) {
            errors <- c(errors,
                "use_idtaxa is TRUE but training_set is NULL"
            )
        } else {
            ts_path <- file.path(config$path, "Data", "Classifiers",
                                 config$training_set)
            if (!file.exists(ts_path)) {
                errors <- c(errors, paste0(
                    "IDTAXA training set not found: ", ts_path
                ))
            }
        }
    }

    if (config$use_blast) {
        if (is.null(config$blast_db)) {
            errors <- c(errors,
                "use_blast is TRUE but blast_db is NULL"
            )
        } else {
            db_dir <- file.path(config$path, "Data", "BlastDBs",
                                config$blast_db)
            if (!dir.exists(db_dir)) {
                errors <- c(errors, paste0(
                    "BLAST database directory not found: ", db_dir
                ))
            }
        }
    }

    # --- External tools ---
    required_tools <- character()
    if (config$use_cutadapt)          required_tools <- c(required_tools, "cutadapt")
    if (config$linkage == "complete")  required_tools <- c(required_tools, "vsearch")
    if (config$linkage == "single")    required_tools <- c(required_tools, "swarm")
    if (config$use_blast)              required_tools <- c(required_tools, "blastn")
    # vsearch is also needed for LULU matchlists
    if (!"vsearch" %in% required_tools) {
        required_tools <- c(required_tools, "vsearch")
    }

    missing_tools <- character()
    for (tool in required_tools) {
        if (Sys.which(tool) == "") {
            missing_tools <- c(missing_tools, tool)
        }
    }
    if (length(missing_tools) > 0) {
        errors <- c(errors, paste0(
            "Required external tools not found on PATH: ",
            paste(missing_tools, collapse = ", "),
            ". Is the conda environment activated?"
        ))
    }

    # --- Report ---
    if (length(warnings) > 0) {
        for (w in warnings) {
            cli::cli_warn(w)
        }
    }

    if (length(errors) > 0) {
        cli::cli_h1("Validation failed")
        for (e in errors) {
            cli::cli_alert_danger(e)
        }
        cli::cli_abort("Fix the above errors before running the pipeline.")
    }

    cli::cli_alert_success("Configuration validated successfully")
    invisible(TRUE)
}
