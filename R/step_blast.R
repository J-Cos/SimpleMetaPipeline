#' Run BLAST taxonomic assignment
#'
#' Performs BLAST searches against a local database and assigns taxonomy
#' at rank-specific similarity thresholds.
#'
#' @param config An `smp_config` object.
#' @param idtaxa_result Result from [step_idtaxa()].
#' @return A list with `SeqDataTable` and `FullBlastOutput`, or `NULL`.
#' @keywords internal
step_blast <- function(config, idtaxa_result) {
    if (!config$use_blast) {
        cli::cli_alert_info("BLAST skipped (use_blast = FALSE)")
        result <- NULL
        smp_cache(result, config, "blast")
        return(result)
    }

    cli::cli_h2("Step 7: BLAST")

    input_table <- idtaxa_result$SeqDataTable
    clustering <- "ESV"

    # Query fasta file
    query_fasta <- paste0(config$dataname, "_", clustering, "_sequences.fasta")
    query_path <- file.path(config$path, "IntermediateOutputs", query_fasta)

    db_path <- file.path(
        config$path, "Data", "BlastDBs",
        config$blast_db, config$blast_db
    )
    output_path <- file.path(
        config$path, "IntermediateOutputs",
        paste0("BlastOutput_", config$blast_db, ".out")
    )

    # Run blastn
    smp_system(
        "blastn",
        c("-db", db_path,
          "-query", query_path,
          "-out", output_path,
          "-outfmt", "'6 qseqid sseqid pident evalue qcovs'",
          "-evalue", "1e-10"),
        step_name = "BLAST"
    )

    # Read output
    blast_output <- utils::read.table(
        output_path, header = FALSE, sep = "\t"
    )
    names(blast_output) <- c(clustering, "SeqID", "Blast_percentIdentical",
                             "Blast_evalue", "Blast_query_coverage")

    # Standardise SeqID casing
    blast_output$SeqID <- stringr::str_to_title(blast_output$SeqID)

    # Load taxonomy mapping
    taxid_path <- file.path(
        config$path, "Data", "BlastDBs",
        paste0(config$blast_db, "_TaxidMap.txt")
    )
    taxid_lines <- readLines(taxid_path)
    taxid_components <- stringr::str_split(taxid_lines, " ", 2)

    taxid_df <- data.frame(
        SeqID = unlist(lapply(taxid_components, `[`, 1)),
        stringsAsFactors = FALSE
    )
    taxa <- unlist(lapply(taxid_components, `[`, 2))
    taxa_ranks <- stringr::str_split(taxa, ";", simplify = TRUE)
    taxa_df <- as.data.frame(taxa_ranks, stringsAsFactors = FALSE)
    names(taxa_df) <- config$blast_desired_ranks
    taxid_df <- cbind(taxid_df, taxa_df)

    # Combine taxonomy with BLAST output
    full_blast <- dplyr::left_join(blast_output, taxid_df, by = "SeqID")
    full_blast <- full_blast[!names(full_blast) == "SeqID"]

    # BUG FIX: Explicitly sort by pident desc, evalue asc before dedup
    full_blast <- full_blast[
        order(-full_blast$Blast_percentIdentical, full_blast$Blast_evalue), ]

    # Pick top hits
    top_hits <- full_blast[!duplicated(full_blast[[clustering]]), ]

    # Apply rank-specific similarity thresholds
    thresholds <- config$blast_thresholds
    for (row_idx in seq_len(nrow(top_hits))) {
        item <- top_hits[row_idx, ]
        # Fill empty assignments
        empty_cells <- which(item == "")
        if (length(empty_cells) > 0) {
            last_non_empty <- max(which(item != ""))
            item[empty_cells] <- paste0("Unknown_", item[last_non_empty])
        }

        # Truncate assignments below threshold
        rank_start <- which(names(item) == config$blast_desired_ranks[1])
        rank_end <- which(names(item) == utils::tail(config$blast_desired_ranks, 1))

        if (length(rank_start) > 0 && length(rank_end) > 0) {
            positive <- item$Blast_percentIdentical > thresholds
            rank_cols <- rank_start:rank_end
            n_positive <- sum(positive[seq_along(rank_cols)])
            if (n_positive < length(rank_cols)) {
                last_positive_rank <- if (n_positive > 0) {
                    item[[rank_cols[n_positive]]]
                } else {
                    "unclassified"
                }
                item[rank_cols[(n_positive + 1):length(rank_cols)]] <-
                    paste0("unclassified_", last_positive_rank)
            }
        }
        top_hits[row_idx, ] <- item
    }

    # Merge with SeqDataTable
    seq_data_table <- merge(
        input_table, top_hits, by = clustering, all = TRUE
    )

    result <- list(
        SeqDataTable = seq_data_table,
        FullBlastOutput = full_blast
    )

    smp_cache(result, config, "blast")
    cli::cli_alert_success("BLAST complete")
    result
}
