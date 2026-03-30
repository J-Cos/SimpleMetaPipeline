#' Run IDTAXA taxonomic classification
#'
#' Classifies sequences using the IDTAXA algorithm from the DECIPHER package.
#' Supports classification of ESVs, OTUs, or curated OTUs.
#'
#' @param config An `smp_config` object.
#' @param input_table The current SeqDataTable.
#' @return A list with `SeqDataTable` (updated with assignments) and
#'   `PlotData` (IDTAXA result object for plotting).
#' @keywords internal
step_idtaxa <- function(config, input_table) {
    if (!config$use_idtaxa) {
        cli::cli_alert_info("IDTAXA skipped (use_idtaxa = FALSE)")
        result <- list(SeqDataTable = input_table, PlotData = NULL)
        smp_cache(result, config, "idtaxa")
        return(result)
    }

    cli::cli_h2("Step 6: IDTAXA")

    # Select sequences to assign
    if (config$seqs_to_assign == "ESVs") {
        dna <- Biostrings::DNAStringSet(input_table[, 2])
        names(dna) <- input_table$ESV
    } else if (config$seqs_to_assign == "OTUs") {
        rows <- input_table$OTUrepresentativeSequence == TRUE
        dna <- Biostrings::DNAStringSet(input_table[rows, 2])
        names(dna) <- input_table$OTU[rows]
    } else if (config$seqs_to_assign == "cOTUs") {
        rows <- input_table$CuratedOTURepresentativeSequence == TRUE
        dna <- Biostrings::DNAStringSet(input_table[rows, 2])
        names(dna) <- input_table$curatedOTU[rows]
    }

    # Load training set
    ts_path <- file.path(config$path, "Data", "Classifiers", config$training_set)
    load(ts_path)  # loads 'trainingSet' into scope

    # Split into chunks
    chunk_size <- 500L
    dna_list <- split(dna, ceiling(seq_along(dna) / chunk_size))
    cli::cli_alert_info(
        "Query sequences split into {length(dna_list)} chunks of {chunk_size}"
    )

    cores <- resolve_cores(config$multithread)

    if (config$idtaxa_parallel) {
        cli::cli_alert("Starting parallel classification with {cores} cores")
        ids <- parallel::mclapply(
            dna_list,
            mc.cores = cores,
            function(chunk) {
                DECIPHER::IdTaxa(
                    chunk, trainingSet,
                    type = "extended",
                    strand = "both",
                    threshold = config$idtaxa_threshold,
                    processors = 1
                )
            }
        )
    } else {
        # Sequential with per-chunk caching for resumability
        ids_cache_path <- file.path(
            config$path, "IntermediateOutputs",
            paste0(config$dataname, "_ids.RDS")
        )
        if (file.exists(ids_cache_path)) {
            ids <- readRDS(ids_cache_path)
        } else {
            ids <- list()
        }

        cli::cli_alert("Starting classification at chunk {length(ids) + 1}")

        for (chunk_num in (length(ids) + 1):length(dna_list)) {
            ids[[chunk_num]] <- DECIPHER::IdTaxa(
                dna_list[[chunk_num]],
                trainingSet,
                type = "extended",
                strand = "both",
                threshold = config$idtaxa_threshold,
                processors = NULL
            )
            cli::cli_alert_success("Chunk {chunk_num} complete")
            saveRDS(ids, ids_cache_path)
        }
    }

    cli::cli_alert_success("All chunks complete")
    ids <- purrr::reduce(ids, c)
    cli::cli_alert_success("All chunks merged")

    # Format results
    idtaxa_df <- format_idtaxa_results(ids, config$desired_ranks)

    # Merge with SeqDataTable
    seq_data_table <- merge(input_table, idtaxa_df, by = "ESV", all = TRUE)

    result <- list(SeqDataTable = seq_data_table, PlotData = ids)
    smp_cache(result, config, "idtaxa")
    cli::cli_alert_success("IDTAXA complete")
    result
}


#' Format IDTAXA results into a data frame
#'
#' Handles both ranked and unranked training sets.
#'
#' @param ids IDTAXA result object.
#' @param desired_ranks Character vector of desired rank names, or `NULL`.
#' @return A data frame with ESV, rank columns, and confidence columns.
#' @keywords internal
format_idtaxa_results <- function(ids, desired_ranks) {
    # Check if rank information is present
    first_item <- as.data.frame(ids[1][[1]])

    if ("rank" %in% colnames(first_item)) {
        # Ranked training set
        reformat_id <- function(id) {
            x <- c(id$taxon, id$confidence)
            names(x) <- c(id$rank, paste0(id$rank, "_confidence"))
            x <- x[c(desired_ranks, paste0(desired_ranks, "_confidence"))]
            names(x) <- c(desired_ranks, paste0(desired_ranks, "_confidence"))
            return(x)
        }

        output_list <- lapply(ids, reformat_id)
        idtaxa_df <- dplyr::bind_rows(lapply(output_list, function(x) {
            as.data.frame(t(x), stringsAsFactors = FALSE)
        }))
        idtaxa_df$ESV <- names(ids)

    } else {
        # Unranked training set
        n_ranks <- max(unlist(lapply(lapply(ids, `[[`, 1), length)))

        extract_from_ids <- function(list_item, category) {
            vals <- list_item[[category]]
            if (length(vals) < n_ranks) {
                vals[(length(vals) + 1):n_ranks] <- NA
            }
            return(vals)
        }

        taxon_list <- lapply(ids, extract_from_ids, category = "taxon")
        taxon_df <- dplyr::bind_rows(lapply(taxon_list, function(x) {
            as.data.frame(t(x), stringsAsFactors = FALSE)
        }))

        conf_list <- lapply(ids, extract_from_ids, category = "confidence")
        conf_df <- dplyr::bind_rows(lapply(conf_list, function(x) {
            as.data.frame(t(x), stringsAsFactors = FALSE)
        }))

        idtaxa_df <- cbind(
            data.frame(ESV = names(ids), stringsAsFactors = FALSE),
            taxon_df, conf_df
        )
        names(idtaxa_df) <- c(
            "ESV",
            paste0("Rank_", seq_len(n_ranks)),
            paste0("Rank_", seq_len(n_ranks), "_Confidence")
        )
    }

    idtaxa_df
}
