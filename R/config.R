#' Create a SimpleMetaPipeline configuration
#'
#' Constructs and validates a configuration object that parameterises the
#' entire pipeline. This replaces the old pattern of defining global variables
#' in a control script.
#'
#' @param path Path to the BioinformaticPipeline_Env working directory.
#' @param dataname Name for this dataset (used to label all outputs).
#' @param data_path Optional path to FASTQs if stored externally. Defaults to
#'   `path`.
#' @param multithread Integer number of threads, or `TRUE` for auto-detect.
#' @param seed Integer seed for reproducibility.
#' @param use_cutadapt Logical; whether to run Cutadapt primer trimming.
#' @param fwd_primer Forward primer sequence.
#' @param rev_primer Reverse primer sequence.
#' @param n_runs Number of sequencing runs.
#' @param trunc_len List of integer vectors; truncation lengths per run.
#' @param trim_left List of integer vectors; bases to trim from start per run.
#' @param max_n Maximum Ns allowed after truncation.
#' @param max_ee Maximum expected errors (forward, reverse).
#' @param trunc_q Truncate reads at first quality score <= this value.
#' @param desired_seq_length_range Optional integer range of sequence lengths to
#'   keep (e.g. `360:270`).
#' @param pool DADA2 pooling: `TRUE`, `FALSE`, or `"pseudo"`.
#' @param mixed_orientation Logical; whether reads are in mixed orientations.
#' @param read_type One of `"paired-end"`, `"single-read"`, or `"premerged"`.
#' @param lulu1_match_rate Percent match rate for LULU stage 1.
#' @param lulu1_min_rel_cooccurrence Minimum relative co-occurrence for LULU
#'   stage 1.
#' @param lulu1_ratio_type Ratio type for LULU stage 1: `"min"` or `"avg"`.
#' @param linkage Clustering linkage: `"complete"` (vsearch) or `"single"`
#'   (swarm).
#' @param differences Number of base differences for swarm clustering.
#' @param similarity_threshold Similarity threshold for vsearch clustering.
#' @param threads Number of threads for clustering tools.
#' @param lulu2_match_rate Percent match rate for LULU stage 2.
#' @param lulu2_min_rel_cooccurrence Minimum relative co-occurrence for LULU
#'   stage 2.
#' @param lulu2_ratio_type Ratio type for LULU stage 2: `"min"` or `"avg"`.
#' @param use_idtaxa Logical; whether to run IDTAXA classification.
#' @param training_set Filename of the IDTAXA training set in
#'   `Data/Classifiers/`.
#' @param desired_ranks Character vector of taxonomic ranks to report.
#' @param seqs_to_assign Which sequences to classify: `"ESVs"`, `"OTUs"`, or
#'   `"cOTUs"`.
#' @param idtaxa_threshold Confidence threshold for IDTAXA assignment (0-100).
#' @param idtaxa_parallel Logical; whether to parallelise IDTAXA classification.
#' @param use_blast Logical; whether to run BLAST.
#' @param blast_db Name of the BLAST database directory.
#' @param blast_thresholds Similarity thresholds per rank for BLAST assignment.
#' @param blast_desired_ranks Character vector of BLAST taxonomic ranks.
#'
#' @return An object of class `smp_config`.
#' @export
smp_config <- function(
    path,
    dataname,
    data_path = NULL,
    multithread = 4L,
    seed = 42L,
    # Cutadapt
    use_cutadapt = FALSE,
    fwd_primer = "ACTG",
    rev_primer = "ACTG",
    # DADA2
    n_runs = 1L,
    trunc_len = list(c(150, 150)),
    trim_left = list(c(0, 0)),
    max_n = 0L,
    max_ee = c(2, 2),
    trunc_q = 2L,
    desired_seq_length_range = NULL,
    pool = "pseudo",
    mixed_orientation = FALSE,
    read_type = c("paired-end", "single-read", "premerged"),
    # LULU 1
    lulu1_match_rate = 90,
    lulu1_min_rel_cooccurrence = 0.95,
    lulu1_ratio_type = c("min", "avg"),
    # Clustering
    linkage = c("complete", "single"),
    differences = 1L,
    similarity_threshold = 0.97,
    threads = 4L,
    # LULU 2
    lulu2_match_rate = 90,
    lulu2_min_rel_cooccurrence = 0.95,
    lulu2_ratio_type = c("min", "avg"),
    # IDTAXA
    use_idtaxa = TRUE,
    training_set = NULL,
    desired_ranks = c("rootrank", "domain", "phylum", "class", "order",
                      "family", "genus"),
    seqs_to_assign = c("ESVs", "OTUs", "cOTUs"),
    idtaxa_threshold = 40,
    idtaxa_parallel = FALSE,
    # BLAST
    use_blast = FALSE,
    blast_db = NULL,
    blast_thresholds = c(0, 0, 80, 85, 90, 95, 97),
    blast_desired_ranks = c("Root", "Kingdom", "Phylum", "Class", "Order",
                            "Family", "Genus", "Species")
) {
    # Match constrained arguments
    read_type <- match.arg(read_type)
    linkage <- match.arg(linkage)
    lulu1_ratio_type <- match.arg(lulu1_ratio_type)
    lulu2_ratio_type <- match.arg(lulu2_ratio_type)
    seqs_to_assign <- match.arg(seqs_to_assign)

    # Resolve data_path default
    if (is.null(data_path)) data_path <- path

    # Resolve multithread
    if (isTRUE(multithread)) {
        multithread <- parallel::detectCores(logical = FALSE)
        cli::cli_alert_info("Auto-detected {multithread} cores")
    }
    multithread <- as.integer(multithread)

    # Build config list
    config <- list(
        path            = normalizePath(path, mustWork = FALSE),
        dataname        = dataname,
        data_path       = normalizePath(data_path, mustWork = FALSE),
        multithread     = multithread,
        seed            = as.integer(seed),
        # Cutadapt
        use_cutadapt    = use_cutadapt,
        fwd_primer      = fwd_primer,
        rev_primer      = rev_primer,
        # DADA2
        n_runs          = as.integer(n_runs),
        trunc_len       = trunc_len,
        trim_left       = trim_left,
        max_n           = as.integer(max_n),
        max_ee          = max_ee,
        trunc_q         = as.integer(trunc_q),
        desired_seq_length_range = desired_seq_length_range,
        pool            = pool,
        mixed_orientation = mixed_orientation,
        read_type       = read_type,
        # LULU 1
        lulu1_match_rate           = lulu1_match_rate,
        lulu1_min_rel_cooccurrence = lulu1_min_rel_cooccurrence,
        lulu1_ratio_type           = lulu1_ratio_type,
        # Clustering
        linkage              = linkage,
        differences          = as.integer(differences),
        similarity_threshold = similarity_threshold,
        threads              = as.integer(threads),
        # LULU 2
        lulu2_match_rate           = lulu2_match_rate,
        lulu2_min_rel_cooccurrence = lulu2_min_rel_cooccurrence,
        lulu2_ratio_type           = lulu2_ratio_type,
        # IDTAXA
        use_idtaxa       = use_idtaxa,
        training_set     = training_set,
        desired_ranks    = desired_ranks,
        seqs_to_assign   = seqs_to_assign,
        idtaxa_threshold = idtaxa_threshold,
        idtaxa_parallel  = idtaxa_parallel,
        # BLAST
        use_blast           = use_blast,
        blast_db            = blast_db,
        blast_thresholds    = blast_thresholds,
        blast_desired_ranks = blast_desired_ranks
    )

    structure(config, class = "smp_config")
}


#' Print method for smp_config
#' @param x An object of class `smp_config`.
#' @param ... Ignored.
#' @export
print.smp_config <- function(x, ...) {
    cli::cli_h1("SimpleMetaPipeline Configuration")
    cli::cli_text("Dataset: {.val {x$dataname}}")
    cli::cli_text("Path: {.path {x$path}}")
    cli::cli_text("Threads: {x$multithread}")
    cli::cli_text("Seed: {x$seed}")
    cli::cli_rule()
    steps <- c()
    if (x$use_cutadapt) steps <- c(steps, "Cutadapt")
    steps <- c(steps, "DADA2", paste0("LULU (ESV, match=", x$lulu1_match_rate, "%)"))
    steps <- c(steps, paste0("Clustering (", x$linkage, ")"))
    steps <- c(steps, paste0("LULU (OTU, match=", x$lulu2_match_rate, "%)"))
    if (x$use_idtaxa) steps <- c(steps, "IDTAXA")
    if (x$use_blast) steps <- c(steps, "BLAST")
    cli::cli_text("Pipeline steps: {.val {paste(steps, collapse = ' \u2192 ')}}")
    invisible(x)
}
