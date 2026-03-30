#' Get all orientations of a primer sequence
#'
#' @param primer A character string representing a DNA primer sequence.
#' @return A named character vector with Forward, Complement, Reverse, and
#'   RevComp orientations.
#' @keywords internal
allOrients <- function(primer) {
    dna <- Biostrings::DNAString(primer)
    orients <- c(
        Forward    = as.character(dna),
        Complement = as.character(Biostrings::complement(dna)),
        Reverse    = as.character(Biostrings::reverse(dna)),
        RevComp    = as.character(Biostrings::reverseComplement(dna))
    )
    return(orients)
}


#' Count primer hits in a FASTQ file
#'
#' @param primer Character string of the primer sequence.
#' @param fn Path to a FASTQ file.
#' @return Integer count of primer hits.
#' @keywords internal
primerHits <- function(primer, fn) {
    nhits <- Biostrings::vcountPattern(
        primer,
        ShortRead::sread(ShortRead::readFastq(fn)),
        fixed = FALSE
    )
    return(sum(nhits > 0))
}


#' Get number of unique sequences from a dada object
#'
#' @param x A dada-class object.
#' @return Integer count.
#' @keywords internal
getN <- function(x) sum(dada2::getUniques(x))


#' Create output file paths for filtered sequences
#'
#' @param config An `smp_config` object.
#' @param sample_names Character vector of sample names.
#' @param suffix File suffix (e.g. `"_F_filt.fastq.gz"`).
#' @param dir_prefix Directory prefix (e.g. `"_filteredsequences"`).
#' @return Character vector of output file paths.
#' @keywords internal
createOutputFilePaths <- function(config, sample_names, suffix, dir_prefix) {
    out_dir <- file.path(
        config$path, "IntermediateOutputs",
        paste0(config$dataname, dir_prefix)
    )
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    file.path(out_dir, paste0(sample_names, suffix))
}


#' Write a clustering summary table
#'
#' @param final_output The final SeqDataTable after LULU stage 2.
#' @param config An `smp_config` object.
#' @keywords internal
writeClusteringTable <- function(final_output, config) {
    n_esvs <- length(unique(final_output$ESV))
    n_otus <- if ("OTU" %in% names(final_output)) {
        length(unique(final_output$OTU))
    } else {
        NA
    }
    n_curated_otus <- if ("curatedOTU" %in% names(final_output)) {
        length(unique(final_output$curatedOTU))
    } else {
        NA
    }

    clustering_table <- data.frame(
        Metric = c("ESVs", "OTUs", "Curated OTUs"),
        Count  = c(n_esvs, n_otus, n_curated_otus)
    )

    utils::write.csv(
        clustering_table,
        file = file.path(
            config$path, "Results",
            paste0(config$dataname, "_ClusteringTable.csv")
        ),
        row.names = FALSE
    )
}


#' Create a FASTA file with abundance annotations for Swarm
#'
#' @param seq_data_table The current SeqDataTable.
#' @param clustering Which clustering level to use (`"ESV"`).
#' @param config An `smp_config` object.
#' @keywords internal
createFastaWithAbundances <- function(seq_data_table, clustering, config) {
    sample_cols <- grep("^Sample_", names(seq_data_table), value = TRUE)
    abundances <- rowSums(seq_data_table[, sample_cols, drop = FALSE])

    seqinr::write.fasta(
        sequences = as.list(paste0(seq_data_table$Sequence, "_", abundances)),
        names     = paste0(seq_data_table$Sequence, "_", abundances),
        file.out  = file.path(
            config$path, "IntermediateOutputs",
            paste0(config$dataname, "_ESVSequencesWithAbundances.fasta")
        )
    )
}


#' Resolve multithread to integer core count
#'
#' @param multithread Integer or TRUE.
#' @return Integer core count.
#' @keywords internal
resolve_cores <- function(multithread) {
    if (isTRUE(multithread)) {
        return(parallel::detectCores(logical = FALSE))
    }
    as.integer(multithread)
}
