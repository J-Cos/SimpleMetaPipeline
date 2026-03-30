#' Run DADA2 denoising
#'
#' Performs quality filtering, error learning, denoising, merging (if
#' paired-end), and chimera removal using DADA2. Handles paired-end,
#' single-read, and premerged read types.
#'
#' @param config An `smp_config` object.
#' @param cutadapt_result Result from [step_cutadapt()].
#' @return A list with `SeqDataTable` and `SecondaryOutputs`.
#' @keywords internal
step_dada2 <- function(config, cutadapt_result) {
    cli::cli_h2("Step 2: DADA2")

    run_names <- list.files(
        file.path(config$data_path, "FASTQs", config$dataname)
    )

    # Check for partially completed runs
    run_esv_tables_path <- file.path(
        config$path, "IntermediateOutputs",
        paste0(config$dataname, "_RunESVtables.RDS")
    )

    if (file.exists(run_esv_tables_path)) {
        RunESVtables <- readRDS(run_esv_tables_path)
        RunDadaPlots <- readRDS(file.path(
            config$path, "IntermediateOutputs",
            paste0(config$dataname, "_RunDadaPlots.RDS")
        ))
        RunDadaTables <- readRDS(file.path(
            config$path, "IntermediateOutputs",
            paste0(config$dataname, "_RunDadaTables.RDS")
        ))
    } else {
        RunESVtables <- list()
        RunDadaPlots <- list()
        RunDadaTables <- list()
    }

    n_completed <- length(RunESVtables)

    if (n_completed < config$n_runs) {
        for (run_idx in (n_completed + 1):config$n_runs) {
            cli::cli_alert_info("Processing Run {run_idx}")
            run_dir <- file.path(
                config$data_path, "FASTQs", config$dataname,
                run_names[run_idx]
            )

            run_result <- dada2_process_run(
                config, run_idx, run_dir, cutadapt_result
            )

            RunESVtables[[run_idx]] <- run_result$esv_table
            RunDadaPlots[[run_idx]] <- run_result$plots
            RunDadaTables[[paste0("Run", run_idx)]] <- run_result$track

            # Cache per-run progress
            saveRDS(RunESVtables, file.path(
                config$path, "IntermediateOutputs",
                paste0(config$dataname, "_RunESVtables.RDS")
            ))
            saveRDS(RunDadaPlots, file.path(
                config$path, "IntermediateOutputs",
                paste0(config$dataname, "_RunDadaPlots.RDS")
            ))
            saveRDS(RunDadaTables, file.path(
                config$path, "IntermediateOutputs",
                paste0(config$dataname, "_RunDadaTables.RDS")
            ))
        }
    }

    # Merge runs
    if (config$n_runs > 1) {
        ESVtable <- dada2::mergeSequenceTables(
            tables = RunESVtables, repeats = "sum"
        )
        DadaPlots <- unlist(RunDadaPlots, recursive = FALSE)
        DadaTables <- do.call("rbind", RunDadaTables)
    } else {
        ESVtable <- RunESVtables[[1]]
        DadaPlots <- RunDadaPlots[[1]]
        DadaTables <- RunDadaTables[[1]]
    }

    SeqLengthDist <- table(nchar(dada2::getSequences(ESVtable)))

    # Chimera removal
    cli::cli_alert("Removing chimeras...")
    ESVtable <- dada2::removeBimeraDenovo(
        ESVtable, method = "consensus",
        multithread = config$multithread, verbose = TRUE
    )

    # Build chimera tracking column
    chimeric_data <- data.frame(
        rowname_for_merging = rownames(ESVtable),
        nonchimeric = rowSums(ESVtable),
        stringsAsFactors = FALSE
    )

    DadaTables <- DadaTables |>
        as.data.frame() |>
        tibble::rownames_to_column("rowname") |>
        tibble::as_tibble() |>
        dplyr::mutate(
            rowname_for_merging = stringr::str_replace(.data$rowname, "_FO", "")
        ) |>
        dplyr::mutate(
            rowname_for_merging = stringr::str_replace(.data$rowname_for_merging, "_RO", "")
        ) |>
        dplyr::left_join(chimeric_data, by = "rowname_for_merging") |>
        dplyr::select(!.data$rowname_for_merging) |>
        dplyr::arrange(.data$rowname)

    # Format ESV table
    ESVtable <- t(ESVtable)  # samples as columns
    ESVtable <- cbind(
        paste0("ESV_", seq_len(nrow(ESVtable))),
        tibble::rownames_to_column(as.data.frame(ESVtable))
    )
    names(ESVtable)[1:2] <- c("ESV", "Sequence")
    names(ESVtable)[3:length(names(ESVtable))] <- paste0(
        "Sample_", names(ESVtable)[3:length(names(ESVtable))]
    )

    result <- list(
        SeqDataTable = ESVtable,
        SecondaryOutputs = list(
            DadaPlots = DadaPlots,
            DadaTables = DadaTables,
            SeqLengthDist = SeqLengthDist
        )
    )

    smp_cache(result, config, "dada2")
    cli::cli_alert_success("DADA2 complete: {nrow(ESVtable)} ESVs")
    result
}


#' Process a single DADA2 run
#'
#' Internal function handling filtering, error learning, denoising, and
#' merging for one sequencing run.
#'
#' @param config An `smp_config` object.
#' @param run_idx Integer run index.
#' @param run_dir Path to the run directory.
#' @param cutadapt_result Result from [step_cutadapt()].
#' @return A list with `esv_table`, `plots`, and `track`.
#' @keywords internal
dada2_process_run <- function(config, run_idx, run_dir, cutadapt_result) {
    if (config$read_type == "paired-end") {
        return(dada2_paired_end(config, run_idx, run_dir, cutadapt_result))
    } else if (config$read_type == "single-read") {
        return(dada2_single_read(config, run_idx, run_dir))
    } else if (config$read_type == "premerged") {
        return(dada2_premerged(config, run_idx, run_dir))
    }
}


#' DADA2 paired-end processing
#' @keywords internal
#' @noRd
dada2_paired_end <- function(config, run_idx, run_dir, cutadapt_result) {
    # Get input files
    if (!config$use_cutadapt || isTRUE(cutadapt_result$skipped)) {
        fnFs <- sort(list.files(run_dir, pattern = "_R1_001.fastq",
                                full.names = TRUE))
        fnRs <- sort(list.files(run_dir, pattern = "_R2_001.fastq",
                                full.names = TRUE))
    } else {
        cut_dir <- file.path(
            config$path, "IntermediateOutputs",
            paste0(config$dataname, "_CutadaptedSeqs")
        )
        fnFs <- sort(list.files(cut_dir, pattern = "_R1_001.fastq",
                                full.names = TRUE))
        fnRs <- sort(list.files(cut_dir, pattern = "_R2_001.fastq",
                                full.names = TRUE))
    }

    sample_names <- paste0(
        sapply(strsplit(basename(fnFs), "_R1_001.fastq"), `[`, 1),
        "__Run", run_idx
    )

    filtFs <- createOutputFilePaths(
        config, sample_names, "_F_filt.fastq.gz", "_filteredsequences"
    )
    filtRs <- createOutputFilePaths(
        config, sample_names, "_R_filt.fastq.gz", "_filteredsequences"
    )

    # Filter and trim
    out <- dada2::filterAndTrim(
        fnFs, filtFs, fnRs, filtRs,
        truncLen = config$trunc_len[[run_idx]],
        trimLeft = config$trim_left[[run_idx]],
        maxN = config$max_n,
        maxEE = config$max_ee,
        truncQ = config$trunc_q,
        rm.phix = TRUE, compress = TRUE,
        multithread = config$multithread, verbose = TRUE
    )
    cli::cli_alert_success("filterAndTrim complete")

    # Learn errors
    errF <- dada2::learnErrors(filtFs, multithread = config$multithread)
    cli::cli_alert_success("learnErrors (forward) complete")
    errR <- dada2::learnErrors(filtRs, multithread = config$multithread)
    cli::cli_alert_success("learnErrors (reverse) complete")

    # Denoise
    dadaFs <- dada2::dada(filtFs, err = errF,
                          multithread = config$multithread, pool = config$pool)
    cli::cli_alert_success("Denoise (forward) complete")
    dadaRs <- dada2::dada(filtRs, err = errR,
                          multithread = config$multithread, pool = config$pool)
    cli::cli_alert_success("Denoise (reverse) complete")

    # Merge
    mergers <- dada2::mergePairs(
        dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE, minOverlap = 5
    )
    cli::cli_alert_success("Merging complete")

    # Diagnostic plots
    plots <- list(
        dada2::plotQualityProfile(fnFs),
        dada2::plotQualityProfile(fnRs),
        dada2::plotQualityProfile(filtFs),
        dada2::plotQualityProfile(filtRs),
        dada2::plotErrors(errF, nominalQ = TRUE),
        dada2::plotErrors(errR, nominalQ = TRUE)
    )

    # Tracking table
    track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN),
                   sapply(mergers, getN))
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                         "merged")
    rownames(track) <- sample_names

    # Build sequence table
    seqtab <- dada2::makeSequenceTable(mergers)
    rownames(seqtab) <- sample_names

    # Length filter
    if (!is.null(config$desired_seq_length_range)) {
        seqtab <- seqtab[, nchar(colnames(seqtab)) %in%
                             config$desired_seq_length_range, drop = FALSE]
    }

    # Handle mixed orientation
    if (config$mixed_orientation) {
        seqtab <- handle_mixed_orientation(seqtab)
    }

    list(esv_table = seqtab, plots = plots, track = track)
}


#' DADA2 single-read processing
#' @keywords internal
#' @noRd
dada2_single_read <- function(config, run_idx, run_dir) {
    fns <- sort(list.files(run_dir, pattern = ".fastq", full.names = TRUE))
    sample_names <- paste0(
        sapply(strsplit(basename(fns), ".fastq"), `[`, 1),
        "__Run", run_idx
    )

    filts <- createOutputFilePaths(
        config, sample_names, "_filt.fastq.gz", "_filteredsequences"
    )

    out <- dada2::filterAndTrim(
        fns, filts,
        truncLen = config$trunc_len[[run_idx]],
        trimLeft = config$trim_left[[run_idx]],
        maxN = config$max_n,
        maxEE = config$max_ee,
        truncQ = config$trunc_q,
        rm.phix = TRUE, compress = TRUE,
        multithread = config$multithread, verbose = TRUE
    )
    cli::cli_alert_success("filterAndTrim complete")

    err <- dada2::learnErrors(filts, multithread = config$multithread)
    cli::cli_alert_success("learnErrors complete")

    # BUG FIX: single-read uses learned errors directly, NOT inflateErr
    dada_result <- dada2::dada(
        filts, err = err,
        multithread = config$multithread, pool = config$pool
    )
    cli::cli_alert_success("Denoise complete")

    plots <- list(
        dada2::plotQualityProfile(fns),
        dada2::plotQualityProfile(filts),
        dada2::plotErrors(err, nominalQ = TRUE)
    )

    track <- cbind(out, sapply(dada_result, getN))
    colnames(track) <- c("input", "filtered", "denoised")
    rownames(track) <- sample_names

    seqtab <- dada2::makeSequenceTable(dada_result)
    rownames(seqtab) <- sample_names

    if (!is.null(config$desired_seq_length_range)) {
        seqtab <- seqtab[, nchar(colnames(seqtab)) %in%
                             config$desired_seq_length_range, drop = FALSE]
    }

    list(esv_table = seqtab, plots = plots, track = track)
}


#' DADA2 premerged read processing
#' @keywords internal
#' @noRd
dada2_premerged <- function(config, run_idx, run_dir) {
    fns <- sort(list.files(run_dir, pattern = ".fastq", full.names = TRUE))
    sample_names <- paste0(
        sapply(strsplit(basename(fns), ".fastq"), `[`, 1),
        "__Run", run_idx
    )

    filts <- createOutputFilePaths(
        config, sample_names, "_filt.fastq.gz", "_filteredsequences"
    )

    out <- dada2::filterAndTrim(
        fns, filts,
        truncLen = config$trunc_len[[run_idx]],
        trimLeft = config$trim_left[[run_idx]],
        maxN = config$max_n,
        maxEE = config$max_ee,
        truncQ = config$trunc_q,
        rm.phix = TRUE, compress = TRUE,
        multithread = config$multithread, verbose = TRUE
    )
    cli::cli_alert_success("filterAndTrim complete")

    err <- dada2::learnErrors(filts, multithread = config$multithread)
    cli::cli_alert_success("learnErrors complete")

    # Premerged reads: inflate errors to correct distortions from premerging
    inflated_err <- dada2::inflateErr(err, 3)

    dada_result <- dada2::dada(
        filts, err = inflated_err,
        multithread = config$multithread, pool = config$pool
    )
    cli::cli_alert_success("Denoise complete")

    plots <- list(
        dada2::plotQualityProfile(fns),
        dada2::plotQualityProfile(filts),
        dada2::plotErrors(inflated_err, nominalQ = TRUE)
    )

    track <- cbind(out, sapply(dada_result, getN))
    colnames(track) <- c("input", "filtered", "denoised")
    rownames(track) <- sample_names

    seqtab <- dada2::makeSequenceTable(dada_result)
    rownames(seqtab) <- sample_names

    if (!is.null(config$desired_seq_length_range)) {
        seqtab <- seqtab[, nchar(colnames(seqtab)) %in%
                             config$desired_seq_length_range, drop = FALSE]
    }

    list(esv_table = seqtab, plots = plots, track = track)
}


#' Handle mixed orientation reads
#'
#' Reverse complements RO (reverse orientation) samples and merges them with
#' FO (forward orientation) samples.
#'
#' @param esv_table A sequence table matrix.
#' @return A merged sequence table with unified orientation.
#' @keywords internal
handle_mixed_orientation <- function(esv_table) {
    # Get reverse-oriented samples
    rev_samples <- esv_table[!grepl("_FO", rownames(esv_table), fixed = TRUE), ,
                             drop = FALSE]
    rev_table <- rev_samples[, colSums(rev_samples) != 0, drop = FALSE]
    reoriented_seqs <- Biostrings::reverseComplement(
        Biostrings::DNAStringSet(colnames(rev_table))
    )
    colnames(rev_table) <- as.character(reoriented_seqs)

    # Get forward-oriented samples
    fwd_samples <- esv_table[!grepl("_RO", rownames(esv_table), fixed = TRUE), ,
                             drop = FALSE]
    fwd_table <- fwd_samples[, colSums(fwd_samples) != 0, drop = FALSE]

    # Clean orientation suffixes
    rownames(rev_table) <- gsub("_RO", "", rownames(rev_table))
    rownames(fwd_table) <- gsub("_FO", "", rownames(fwd_table))

    # Merge
    dada2::mergeSequenceTables(fwd_table, rev_table, repeats = "sum")
}
