#' Run OTU clustering
#'
#' Clusters ESVs into OTUs using either VSEARCH (complete linkage) or Swarm
#' (single linkage), then merges results back into the SeqDataTable.
#'
#' @param config An `smp_config` object.
#' @param dada2_result Result from [step_dada2()].
#' @param lulu1_result Result from [step_lulu()] stage 1.
#' @return The updated SeqDataTable with OTU columns added.
#' @keywords internal
step_cluster <- function(config, dada2_result, lulu1_result) {
    cli::cli_h2("Step 4: Clustering ({config$linkage} linkage)")

    if (config$linkage == "single") {
        createFastaWithAbundances(lulu1_result, "ESV", config)
        result <- run_swarm(config, lulu1_result)
    } else if (config$linkage == "complete") {
        result <- run_vsearch(config, lulu1_result)
    }

    smp_cache(result, config, "cluster")
    n_otus <- length(unique(result$OTU))
    cli::cli_alert_success("Clustering complete: {n_otus} OTUs")
    result
}


#' Run VSEARCH clustering
#' @keywords internal
run_vsearch <- function(config, table_to_merge) {
    cores <- resolve_cores(config$multithread)
    clustering <- "ESV"

    # Get representative sequences
    if ("CuratedESVRepresentativeSequence" %in% names(table_to_merge)) {
        seqs <- table_to_merge$Sequence[
            table_to_merge$CuratedESVRepresentativeSequence
        ]
        clustering <- "curatedESV"
    } else {
        seqs <- table_to_merge$Sequence
    }

    # Write FASTA
    fasta_file <- paste0(config$dataname, "_", clustering,
                         "_sequences_NamedAsSeqs.fasta")
    fasta_path <- file.path(config$path, "IntermediateOutputs", fasta_file)
    seqinr::write.fasta(
        sequences = as.list(seqs),
        names = seqs,
        file.out = fasta_path
    )

    # Run vsearch
    otu_dir <- file.path(
        config$path, "IntermediateOutputs",
        paste0(config$dataname, "_vsearchOTUs")
    )
    centroid_path <- file.path(
        config$path, "IntermediateOutputs",
        paste0(config$dataname, "_centroidsequences.fasta")
    )
    uc_path <- file.path(
        config$path, "IntermediateOutputs",
        paste0(config$dataname, "_vsearchdata")
    )

    if (!dir.exists(otu_dir)) dir.create(otu_dir, recursive = TRUE)

    smp_system(
        "vsearch",
        c("--cluster_size", fasta_path,
          "--centroids", centroid_path,
          "--clusterout_id",
          "--clusters", file.path(otu_dir, "Cluster"),
          "--uc", uc_path,
          "--threads", cores,
          "--id", config$similarity_threshold,
          "--relabel", "OTU"),
        step_name = "VSEARCH clustering"
    )

    # Parse results
    vsearch_data <- utils::read.table(uc_path)
    vsearch_data <- vsearch_data[vsearch_data$V1 != "C", c(1, 2, 9)]

    otu_labeled <- data.frame(
        OTU = paste0("OTU_", vsearch_data$V2),
        Sequence = vsearch_data$V9,
        OTUrepresentativeSequence = (vsearch_data$V1 == "S"),
        stringsAsFactors = FALSE
    )

    merge_clustering_results(table_to_merge, otu_labeled)
}


#' Run Swarm clustering
#' @keywords internal
run_swarm <- function(config, table_to_merge) {
    fasta_path <- file.path(
        config$path, "IntermediateOutputs",
        paste0(config$dataname, "_ESVSequencesWithAbundances.fasta")
    )

    smp_system(
        "swarm",
        c("-d", config$differences,
          "-f",
          "-t", config$threads,
          fasta_path,
          "-s", file.path(config$path, "IntermediateOutputs",
                          paste0(config$dataname, "_SwarmStats.csv")),
          "-w", file.path(config$path, "IntermediateOutputs",
                          paste0(config$dataname, "SwarmSequences.fasta")),
          "-o", file.path(config$path, "IntermediateOutputs",
                          paste0(config$dataname, "_SwarmAllAmplicons.txt")),
          "-i", file.path(config$path, "IntermediateOutputs",
                          paste0(config$dataname, "_SwarmInternalStructures.txt"))),
        step_name = "Swarm clustering"
    )

    # Parse swarm output
    all_amplicons <- readLines(file.path(
        config$path, "IntermediateOutputs",
        paste0(config$dataname, "_SwarmAllAmplicons.txt")
    ))

    seq_data_list <- list()
    for (i in seq_along(all_amplicons)) {
        otu <- data.frame(
            Sequence = unlist(strsplit(all_amplicons[i], " ")),
            OTU = paste0("OTU_", i),
            stringsAsFactors = FALSE
        )
        otu$Sequence <- sapply(strsplit(otu$Sequence, "_"), `[`, 1)
        otu$OTUrepresentativeSequence <- FALSE
        otu$OTUrepresentativeSequence[1] <- TRUE
        seq_data_list[[i]] <- otu
    }
    otu_labeled <- do.call(rbind, seq_data_list)

    merge_clustering_results(table_to_merge, otu_labeled)
}


#' Merge clustering results into SeqDataTable
#'
#' Combines OTU assignments with the existing SeqDataTable, backfilling
#' assignments for ESVs that were removed by LULU curation.
#'
#' @param table_to_merge The current SeqDataTable.
#' @param otu_labeled Data frame with Sequence, OTU, OTUrepresentativeSequence.
#' @return The merged SeqDataTable.
#' @keywords internal
merge_clustering_results <- function(table_to_merge, otu_labeled) {
    result <- merge(table_to_merge, otu_labeled, by = "Sequence", all = TRUE)

    # Backfill OTU assignments for ESVs grouped by LULU into curated ESVs
    removed_esvs <- which(is.na(result[["OTU"]]))
    if (length(removed_esvs) > 0 && "curatedESV" %in% names(result)) {
        for (i in removed_esvs) {
            parent_row <- which(result$ESV == result[i, ]$curatedESV)
            if (length(parent_row) > 0) {
                result[i, "OTU"] <- result[parent_row[1], "OTU"]
                result[i, "OTUrepresentativeSequence"] <- FALSE
            }
        }
    }

    result
}
