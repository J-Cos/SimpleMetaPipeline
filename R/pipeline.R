#' Run the SimpleMetaPipeline
#'
#' Executes the full metabarcoding pipeline: Cutadapt → DADA2 → LULU₁ →
#' Clustering → LULU₂ → IDTAXA → BLAST → Results. Automatically resumes
#' from the last completed step if intermediate outputs exist.
#'
#' @param config An `smp_config` object created by [smp_config()].
#'
#' @return Invisibly returns the cache list containing all step results.
#' @export
#'
#' @examples
#' \dontrun{
#' config <- smp_config(
#'     path = "../BioinformaticPipeline_Env",
#'     dataname = "Example",
#'     n_runs = 2L,
#'     trunc_len = list(c(150, 150), c(150, 150)),
#'     trim_left = list(c(0, 0), c(0, 0)),
#'     training_set = "GTDB_r207-mod_April2022.RData"
#' )
#' smp_run(config)
#' }
smp_run <- function(config) {
    stopifnot(inherits(config, "smp_config"))

    cli::cli_h1("SimpleMetaPipeline v3.0.0")
    cli::cli_text("Dataset: {.val {config$dataname}}")
    cli::cli_text("Read type: {.val {config$read_type}}")
    cli::cli_text("Threads: {config$multithread}")
    cli::cli_rule()

    # Validate before starting
    smp_validate(config)
    set.seed(config$seed)

    # Load any cached results
    cache <- smp_load_cache(config)

    # Step 1: Cutadapt
    if (is.null(cache$cutadapt)) {
        cache$cutadapt <- step_cutadapt(config)
    } else {
        cli::cli_alert_info("Cutadapt: using cached result")
    }

    # Step 2: DADA2
    if (is.null(cache$dada2)) {
        cache$dada2 <- step_dada2(config, cache$cutadapt)
    } else {
        cli::cli_alert_info("DADA2: using cached result")
    }

    # Step 3: LULU stage 1 (ESV curation)
    if (is.null(cache$lulu1)) {
        cache$lulu1 <- step_lulu(
            config, cache$dada2$SeqDataTable, cache$dada2, stage = 1
        )
    } else {
        cli::cli_alert_info("LULU (ESV): using cached result")
    }

    # Step 4: Clustering
    if (is.null(cache$cluster)) {
        cache$cluster <- step_cluster(config, cache$dada2, cache$lulu1)
    } else {
        cli::cli_alert_info("Clustering: using cached result")
    }

    # Step 5: LULU stage 2 (OTU curation)
    if (is.null(cache$lulu2)) {
        cache$lulu2 <- step_lulu(
            config, cache$cluster, cache$dada2, stage = 2
        )
    } else {
        cli::cli_alert_info("LULU (OTU): using cached result")
    }

    # Step 6: IDTAXA
    if (is.null(cache$idtaxa)) {
        cache$idtaxa <- step_idtaxa(config, cache$lulu2)
    } else {
        cli::cli_alert_info("IDTAXA: using cached result")
    }

    # Step 7: BLAST
    if (is.null(cache$blast)) {
        cache$blast <- step_blast(config, cache$idtaxa)
    } else {
        cli::cli_alert_info("BLAST: using cached result")
    }

    # Step 8: Save results
    step_save_results(config, cache)

    cli::cli_rule()
    cli::cli_alert_success("Pipeline complete for {.val {config$dataname}}")

    invisible(cache)
}


#' Save final pipeline results
#'
#' Writes the final SeqDataTable, diagnostic plots, tracking tables, and
#' clustering summaries to the Results directory.
#'
#' @param config An `smp_config` object.
#' @param cache The cache list from [smp_run()].
#' @keywords internal
step_save_results <- function(config, cache) {
    cli::cli_h2("Step 8: Saving results")

    results_dir <- file.path(config$path, "Results")

    # Save DADA2 diagnostic plots (if available and compatible)
    dada_plots <- NULL
    if (is.list(cache$dada2) && "SecondaryOutputs" %in% names(cache$dada2)) {
        dada_plots <- cache$dada2$SecondaryOutputs$DadaPlots
    }
    if (!is.null(dada_plots)) {
        tryCatch({
            glist <- lapply(dada_plots, ggplot2::ggplotGrob)
            ggplot2::ggsave(
                file.path(results_dir,
                           paste0(config$dataname, "_DadaPlots.pdf")),
                gridExtra::marrangeGrob(glist, nrow = 1, ncol = 1)
            )
        }, error = function(e) {
            cli::cli_warn(c(
                "Could not save DADA2 plots (possibly from an older ggplot2 version)",
                "i" = conditionMessage(e)
            ))
        })
    }

    # Save tracking tables (if available)
    if (!is.null(cache$cutadapt) && !isTRUE(cache$cutadapt$skipped)) {
        tryCatch({
            utils::write.csv(
                cache$cutadapt$PrimerCountSamplesList,
                file = file.path(results_dir,
                                 paste0(config$dataname, "_PrimerCounts.csv"))
            )
            utils::write.csv(
                cache$cutadapt$PrimerCountCutadaptedSamplesList,
                file = file.path(results_dir,
                                 paste0(config$dataname,
                                        "_PrimerCountsAfterCutadapt.csv"))
            )
        }, error = function(e) {
            cli::cli_warn("Could not save cutadapt primer counts: {conditionMessage(e)}")
        })
    }

    # DADA2 tracking table
    dada_tables <- NULL
    if (is.list(cache$dada2) && "SecondaryOutputs" %in% names(cache$dada2)) {
        dada_tables <- cache$dada2$SecondaryOutputs$DadaTables
    }
    if (!is.null(dada_tables)) {
        utils::write.csv(
            dada_tables,
            file = file.path(results_dir,
                             paste0(config$dataname, "_DadaTable.csv"))
        )
    }

    # Sequence length distribution
    seq_len_dist <- NULL
    if (is.list(cache$dada2) && "SecondaryOutputs" %in% names(cache$dada2)) {
        seq_len_dist <- cache$dada2$SecondaryOutputs$SeqLengthDist
    }
    if (!is.null(seq_len_dist)) {
        utils::write.csv(
            seq_len_dist,
            file = file.path(results_dir,
                             paste0(config$dataname,
                                    "_DadaSeqLengthDistribution.csv"))
        )
    }

    # Save clustering table
    tryCatch(
        writeClusteringTable(cache$lulu2, config),
        error = function(e) {
            cli::cli_warn("Could not write clustering table: {conditionMessage(e)}")
        }
    )

    # Determine final table — handle both legacy (data.frame) and new (list)
    final_table <- NULL
    if (!is.null(cache$blast)) {
        if (is.data.frame(cache$blast)) {
            final_table <- cache$blast
        } else if ("SeqDataTable" %in% names(cache$blast)) {
            final_table <- cache$blast$SeqDataTable
        }
    }
    if (is.null(final_table) && !is.null(cache$idtaxa)) {
        if (is.data.frame(cache$idtaxa)) {
            final_table <- cache$idtaxa
        } else if ("SeqDataTable" %in% names(cache$idtaxa)) {
            final_table <- cache$idtaxa$SeqDataTable
        }
    }
    if (is.null(final_table)) {
        final_table <- cache$lulu2
    }

    saveRDS(
        final_table,
        file = file.path(results_dir,
                         paste0(config$dataname, "_SeqDataTable.RDS"))
    )
    utils::write.csv(
        final_table,
        file = file.path(results_dir,
                         paste0(config$dataname, "_SeqDataTable.csv"))
    )

    # Save taxonomy assignment plot if available
    if (is.list(cache$idtaxa) && "PlotData" %in% names(cache$idtaxa) &&
        !is.null(cache$idtaxa$PlotData)) {
        tryCatch({
            grDevices::pdf(
                file = file.path(results_dir,
                                 paste0(config$dataname,
                                        "_TaxaAssignment.pdf"))
            )
            plot(cache$idtaxa$PlotData)
            grDevices::dev.off()
        }, error = function(e) {
            try(grDevices::dev.off(), silent = TRUE)
            cli::cli_warn("Could not save taxonomy plot: {conditionMessage(e)}")
        })
    }

    cli::cli_alert_success(
        "Results saved to {.path {results_dir}}"
    )
}


#' Resume a pipeline run
#'
#' Convenience wrapper around [smp_run()] with a clear name indicating
#' that the pipeline will resume from cached intermediate outputs.
#'
#' @param config An `smp_config` object.
#' @return Invisibly returns the cache list.
#' @export
smp_resume <- function(config) {
    cli::cli_alert_info("Resuming pipeline from last completed step")
    smp_run(config)
}
