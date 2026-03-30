#' Run Cutadapt primer trimming
#'
#' Removes primer sequences from paired-end reads using Cutadapt. If
#' `config$use_cutadapt` is `FALSE`, returns a list with `skipped = TRUE`.
#'
#' @param config An `smp_config` object.
#' @return A list with `PrimerCountSamplesList` and
#'   `PrimerCountCutadaptedSamplesList`, or `list(skipped = TRUE)`.
#' @keywords internal
step_cutadapt <- function(config) {
    if (!config$use_cutadapt) {
        cli::cli_alert_info("Cutadapt skipped (use_cutadapt = FALSE)")
        result <- list(
            skipped = TRUE,
            PrimerCountSamplesList = data.frame(),
            PrimerCountCutadaptedSamplesList = data.frame()
        )
        smp_cache(result, config, "cutadapt")
        return(result)
    }

    cli::cli_h2("Step 1: Cutadapt")

    run_names <- list.files(
        file.path(config$data_path, "FASTQs", config$dataname)
    )

    # Get all primer orientations
    FWD_orients <- allOrients(config$fwd_primer)
    REV_orients <- allOrients(config$rev_primer)

    # Process each run
    all_primer_counts <- list()
    all_cutadapted_counts <- list()

    for (run_idx in seq_len(config$n_runs)) {
        run_dir <- file.path(
            config$data_path, "FASTQs", config$dataname, run_names[run_idx]
        )

        fnFs <- sort(list.files(run_dir, pattern = "_R1_001.fastq",
                                full.names = TRUE))
        fnRs <- sort(list.files(run_dir, pattern = "_R2_001.fastq",
                                full.names = TRUE))

        sample_names <- sapply(
            strsplit(basename(fnFs), "_R1_001.fastq"), `[`, 1
        )

        # Pre-filter to remove Ns
        filtFs <- createOutputFilePaths(
            config, sample_names, "_F_filt.fastq.gz", "_prefilteredsequences"
        )
        filtRs <- createOutputFilePaths(
            config, sample_names, "_R_filt.fastq.gz", "_prefilteredsequences"
        )

        dada2::filterAndTrim(
            fnFs, filtFs, fnRs, filtRs,
            maxN = 0, rm.phix = TRUE, compress = TRUE,
            multithread = config$multithread, verbose = TRUE
        )

        # Count primers before cutadapt
        primer_counts <- list()
        for (i in seq_along(filtFs)) {
            primer_counts[[i]] <- rbind(
                FWD.ForwardReads = sapply(FWD_orients, primerHits, fn = filtFs[[i]]),
                FWD.ReverseReads = sapply(FWD_orients, primerHits, fn = filtRs[[i]]),
                REV.ForwardReads = sapply(REV_orients, primerHits, fn = filtFs[[i]]),
                REV.ReverseReads = sapply(REV_orients, primerHits, fn = filtRs[[i]])
            )
        }
        all_primer_counts <- c(all_primer_counts, primer_counts)

        # Run Cutadapt
        path_cut <- file.path(
            config$path, "IntermediateOutputs",
            paste0(config$dataname, "_CutadaptedSeqs")
        )
        if (!dir.exists(path_cut)) dir.create(path_cut)

        fnFs_cut <- file.path(path_cut, basename(fnFs))
        fnRs_cut <- file.path(path_cut, basename(fnRs))

        FWD_RC <- dada2:::rc(config$fwd_primer)
        REV_RC <- dada2:::rc(config$rev_primer)

        R1_flags <- paste("-g", config$fwd_primer, "-a", REV_RC)
        R2_flags <- paste("-G", config$rev_primer, "-A", FWD_RC)

        for (i in seq_along(fnFs)) {
            smp_system(
                "cutadapt",
                c(R1_flags, R2_flags, "-n", "2",
                  "-o", fnFs_cut[i], "-p", fnRs_cut[i],
                  filtFs[i], filtRs[i]),
                step_name = paste0("Cutadapt sample ", i)
            )
        }

        # Count primers after cutadapt
        cutadapted_counts <- list()
        for (i in seq_along(filtFs)) {
            cutadapted_counts[[i]] <- rbind(
                FWD.ForwardReads = sapply(FWD_orients, primerHits, fn = fnFs_cut[[i]]),
                FWD.ReverseReads = sapply(FWD_orients, primerHits, fn = fnRs_cut[[i]]),
                REV.ForwardReads = sapply(REV_orients, primerHits, fn = fnFs_cut[[i]]),
                REV.ReverseReads = sapply(REV_orients, primerHits, fn = fnRs_cut[[i]])
            )
        }
        all_cutadapted_counts <- c(all_cutadapted_counts, cutadapted_counts)
    }

    result <- list(
        skipped = FALSE,
        PrimerCountSamplesList = all_primer_counts,
        PrimerCountCutadaptedSamplesList = all_cutadapted_counts
    )

    smp_cache(result, config, "cutadapt")
    cli::cli_alert_success("Cutadapt complete")
    result
}
