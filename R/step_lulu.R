#' Run LULU curation
#'
#' Performs LULU post-clustering curation to merge co-occurring sequences that
#' are likely to be errors. Can be run at ESV level (stage 1) or OTU level
#' (stage 2).
#'
#' @param config An `smp_config` object.
#' @param input_table The current SeqDataTable to curate.
#' @param dada2_result Result from [step_dada2()] (needed for sample columns).
#' @param stage Integer, 1 (ESV curation) or 2 (OTU curation).
#' @return The updated SeqDataTable with curated columns added.
#' @keywords internal
step_lulu <- function(config, input_table, dada2_result, stage) {
    if (stage == 1) {
        cli::cli_h2("Step 3: LULU (ESV curation)")
        clustering <- "ESV"
        match_rate <- config$lulu1_match_rate
        min_rel_co <- config$lulu1_min_rel_cooccurrence
        ratio_type <- config$lulu1_ratio_type
        cache_name <- "lulu1"
    } else {
        cli::cli_h2("Step 5: LULU (OTU curation)")
        clustering <- "OTU"
        match_rate <- config$lulu2_match_rate
        min_rel_co <- config$lulu2_min_rel_cooccurrence
        ratio_type <- config$lulu2_ratio_type
        cache_name <- "lulu2"
    }

    # Create OTU table for LULU
    otu_table_lulu <- create_otu_table_for_lulu(
        input_table, dada2_result, clustering
    )

    # Create match list
    match_list_lulu <- create_matchlist_for_lulu(
        input_table, config, match_rate, clustering
    )

    # Run LULU
    if (!requireNamespace("lulu", quietly = TRUE)) {
        cli::cli_abort("Package {.pkg lulu} is required. Install with: devtools::install_github('tobiasgf/lulu')")
    }

    curated_result <- lulu::lulu(
        otu_table_lulu, match_list_lulu,
        minimum_ratio_type = ratio_type,
        minimum_match = match_rate,
        minimum_relative_cooccurence = min_rel_co
    )

    # Extract curation mapping
    curated_map <- curated_result$otu_map
    curated_map[[clustering]] <- rownames(curated_map)
    curated_map <- curated_map[, c("parent_id", "curated", clustering)]
    names(curated_map) <- c(
        paste0("curated", clustering),
        paste0("Curated", clustering, "RepresentativeSequence"),
        clustering
    )
    curated_map[[paste0("Curated", clustering, "RepresentativeSequence")]] <-
        curated_map[[paste0("Curated", clustering, "RepresentativeSequence")]] == "parent"

    # Merge with SeqDataTable
    result <- merge(input_table, curated_map, by = clustering, all = TRUE)

    # For OTU stage: ESVs that aren't OTU representatives can't be cOTU reps
    if (clustering == "OTU") {
        result$CuratedOTURepresentativeSequence <-
            result$CuratedOTURepresentativeSequence == TRUE &
            result$OTUrepresentativeSequence == TRUE
    }

    smp_cache(result, config, cache_name)
    cli::cli_alert_success(
        "LULU ({clustering}) complete: {sum(curated_result$otu_map$curated == 'parent')} curated {tolower(clustering)}s"
    )
    result
}


#' Create OTU/ESV table formatted for LULU input
#' @keywords internal
create_otu_table_for_lulu <- function(input_table, dada2_result, clustering) {
    sample_cols <- grep("^Sample_", names(dada2_result$SeqDataTable),
                        value = TRUE)

    if (clustering == "ESV") {
        otu_table <- input_table |>
            dplyr::group_by(.data$ESV) |>
            dplyr::select(dplyr::all_of(sample_cols)) |>
            dplyr::summarise(dplyr::across(dplyr::everything(), sum))
    } else if (clustering == "OTU") {
        otu_table <- input_table |>
            dplyr::group_by(.data$OTU) |>
            dplyr::select(dplyr::all_of(sample_cols)) |>
            dplyr::summarise(dplyr::across(dplyr::everything(), sum))
    }

    otu_table <- as.data.frame(otu_table)
    rownames(otu_table) <- otu_table[[clustering]]
    otu_table <- otu_table[, -1, drop = FALSE]

    otu_table
}


#' Create a match list for LULU using vsearch
#' @keywords internal
create_matchlist_for_lulu <- function(input_table, config, match_rate,
                                      clustering) {
    cores <- resolve_cores(config$multithread)

    if (clustering == "ESV") {
        seqs <- input_table
    } else if (clustering == "OTU") {
        seqs <- input_table[input_table$OTUrepresentativeSequence == TRUE, ]
    }

    # Write FASTA
    fasta_file <- paste0(config$dataname, "_", clustering, "_sequences.fasta")
    fasta_path <- file.path(config$path, "IntermediateOutputs", fasta_file)
    seqinr::write.fasta(
        sequences = as.list(seqs$Sequence),
        names = seqs[[clustering]],
        file.out = fasta_path
    )

    # Run vsearch to create match list
    matchlist_file <- paste0(config$dataname, "_", clustering,
                             "_match_list.txt")
    matchlist_path <- file.path(config$path, "IntermediateOutputs",
                                matchlist_file)

    smp_system(
        "vsearch",
        c("--usearch_global", fasta_path,
          "--db", fasta_path,
          "--threads", cores,
          "--self",
          "--id", paste0(".", match_rate),
          "--iddef", "1",
          "--userout", matchlist_path,
          "-userfields", "query+target+id",
          "--maxaccepts", "0",
          "--maxrejects", "0",
          "--query_cov", ".9",
          "--maxhits", "10"),
        step_name = paste("LULU matchlist", clustering)
    )

    match_list <- utils::read.table(matchlist_path, sep = "\t", header = FALSE)
    match_list
}
