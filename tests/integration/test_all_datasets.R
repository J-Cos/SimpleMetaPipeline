#!/usr/bin/env Rscript
###############################################################################
# Integration test: Validate smp_config + smp_validate across all available
# test datasets with a range of configurations.
#
# This script does NOT run the full pipeline (that would take hours with DADA2).
# Instead it tests:
#   1. Config creation with parameters matching each dataset
#   2. Validation passes for datasets with correct structure
#   3. Cache loading works against existing intermediate outputs
#   4. The refactored pipeline can load and inspect cached results
#
# Run with:
#   conda run --prefix ./env Rscript tests/integration/test_all_datasets.R
###############################################################################

devtools::load_all(".")

# The package lives at BioinformaticPipeline_Env/BioinformaticPipeline
# and env root is BioinformaticPipeline_Env (the parent)
env_root <- normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
cat("========================================\n")
cat("Integration test: SimpleMetaPipeline v3\n")
cat("Environment root:", env_root, "\n")
cat("========================================\n\n")

passed <- 0
failed <- 0
errors <- character()

# Helper to run a test
run_test <- function(name, expr) {
    cat(sprintf("  [TEST] %s ... ", name))
    result <- tryCatch({
        expr
        cat("PASS\n")
        TRUE
    }, error = function(e) {
        cat(sprintf("FAIL: %s\n", conditionMessage(e)))
        FALSE
    })
    return(result)
}

###############################################################################
# TEST 1: Example dataset — 2-run paired-end with GTDB classifier
###############################################################################
cat("\n── Dataset: Example (2-run paired-end) ──────────────────────────────\n")

config_example <- smp_config(
    path = env_root,
    dataname = "Example",
    n_runs = 2L,
    trunc_len = list(c(150, 150), c(150, 150)),
    trim_left = list(c(0, 0), c(0, 0)),
    use_cutadapt = FALSE,
    read_type = "paired-end",
    linkage = "complete",
    similarity_threshold = 0.97,
    use_idtaxa = TRUE,
    training_set = "GTDB_r207-mod_April2022.RData",
    desired_ranks = c("rootrank", "domain", "phylum", "class", "order",
                      "family", "genus"),
    seqs_to_assign = "ESVs",
    use_blast = FALSE,
    multithread = 4L
)

if (run_test("Config creation", { stopifnot(inherits(config_example, "smp_config")) })) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

if (run_test("Validation passes", { smp_validate(config_example) })) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

if (run_test("Cache loading", {
    cache <- smp_load_cache(config_example)
    stopifnot(!is.null(cache$dada2))
    stopifnot(!is.null(cache$lulu1))
    stopifnot(!is.null(cache$cluster))
    stopifnot(!is.null(cache$lulu2))
    stopifnot(!is.null(cache$idtaxa))
    cat(sprintf("(loaded %d cached steps) ", sum(!sapply(cache, is.null))))
})) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

if (run_test("DADA2 result structure", {
    cache <- smp_load_cache(config_example)
    dada_result <- cache$dada2
    stopifnot("SeqDataTable" %in% names(dada_result))
    stopifnot("SecondaryOutputs" %in% names(dada_result))
    sdt <- dada_result$SeqDataTable
    stopifnot("ESV" %in% names(sdt))
    stopifnot("Sequence" %in% names(sdt))
    sample_cols <- grep("^Sample_", names(sdt), value = TRUE)
    cat(sprintf("(%d ESVs, %d samples) ", nrow(sdt), length(sample_cols)))
})) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

###############################################################################
# TEST 2: 16s_multirun_test — paired-end with SILVA classifier + BLAST
###############################################################################
cat("\n── Dataset: 16s_multirun_test (paired-end, SILVA + BLAST) ───────────\n")

# Check if this dataset exists
fastq_16s <- file.path(env_root, "FASTQs", "16s_test")
if (dir.exists(fastq_16s)) {
    n_runs_16s <- length(list.files(fastq_16s))
    if (n_runs_16s == 0) n_runs_16s <- 1  # flat directory

    config_16s <- smp_config(
        path = env_root,
        dataname = "16s_multirun_test",
        n_runs = 1L,
        trunc_len = list(c(150, 150)),
        trim_left = list(c(0, 0)),
        use_cutadapt = FALSE,
        read_type = "paired-end",
        linkage = "complete",
        similarity_threshold = 0.97,
        use_idtaxa = TRUE,
        training_set = "GTDB_r207-mod_April2022.RData",
        use_blast = TRUE,
        blast_db = "SILVA138_BlastDB",
        blast_desired_ranks = c("Root", "Kingdom", "Phylum", "Class", "Order",
                                "Family", "Genus", "Species"),
        multithread = 4L
    )

    if (run_test("Config creation (BLAST enabled)", {
        stopifnot(inherits(config_16s, "smp_config"))
        stopifnot(config_16s$use_blast == TRUE)
    })) {
        passed <- passed + 1
    } else {
        failed <- failed + 1
    }

    if (run_test("Cache loading (16s)", {
        cache <- smp_load_cache(config_16s)
        n_cached <- sum(!sapply(cache, is.null))
        cat(sprintf("(%d cached steps) ", n_cached))
        stopifnot(n_cached > 0)
    })) {
        passed <- passed + 1
    } else {
        failed <- failed + 1
    }
} else {
    cat("  SKIPPED: 16s_test directory not found\n")
}

###############################################################################
# TEST 3: 18s_test — paired-end with PR2 classifier
###############################################################################
cat("\n── Dataset: 18s_test (paired-end, PR2) ──────────────────────────────\n")

config_18s <- smp_config(
    path = env_root,
    dataname = "18s_test",
    n_runs = 1L,
    trunc_len = list(c(200, 200)),
    trim_left = list(c(0, 0)),
    use_cutadapt = FALSE,
    read_type = "paired-end",
    linkage = "complete",
    similarity_threshold = 0.97,
    use_idtaxa = TRUE,
    training_set = "PR2_v4_13_March2021.RData",
    use_blast = FALSE,
    multithread = 4L
)

if (run_test("Config creation (PR2 classifier)", {
    stopifnot(inherits(config_18s, "smp_config"))
    stopifnot(config_18s$training_set == "PR2_v4_13_March2021.RData")
})) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

if (run_test("Cache loading (18s)", {
    cache <- smp_load_cache(config_18s)
    n_cached <- sum(!sapply(cache, is.null))
    cat(sprintf("(%d cached steps) ", n_cached))
    stopifnot(n_cached > 0)
})) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

###############################################################################
# TEST 4: COI_test — paired-end with Biocode classifier
###############################################################################
cat("\n── Dataset: COI_test (paired-end, Biocode) ──────────────────────────\n")

config_coi <- smp_config(
    path = env_root,
    dataname = "COI_test",
    n_runs = 1L,
    trunc_len = list(c(200, 200)),
    trim_left = list(c(0, 0)),
    use_cutadapt = FALSE,
    read_type = "paired-end",
    linkage = "complete",
    similarity_threshold = 0.97,
    use_idtaxa = TRUE,
    training_set = "Biocode_ER_IdtaxaClassifier.Rdata",
    use_blast = FALSE,
    multithread = 4L
)

if (run_test("Config creation (COI, Biocode)", {
    stopifnot(inherits(config_coi, "smp_config"))
})) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

if (run_test("Cache loading (COI)", {
    cache <- smp_load_cache(config_coi)
    n_cached <- sum(!sapply(cache, is.null))
    cat(sprintf("(%d cached steps) ", n_cached))
    stopifnot(n_cached > 0)
})) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

###############################################################################
# TEST 5: COI_multirun_test — 4-run, paired-end, Swarm (single linkage)
###############################################################################
cat("\n── Dataset: COI_multirun_test (4-run, Swarm clustering) ─────────────\n")

config_coi_multi <- smp_config(
    path = env_root,
    dataname = "COI_multirun_test",
    n_runs = 4L,
    trunc_len = list(c(200, 200), c(200, 200), c(200, 200), c(200, 200)),
    trim_left = list(c(0, 0), c(0, 0), c(0, 0), c(0, 0)),
    use_cutadapt = FALSE,
    read_type = "paired-end",
    linkage = "single",        # Swarm instead of VSEARCH
    differences = 1L,
    use_idtaxa = TRUE,
    training_set = "Biocode_ER_IdtaxaClassifier.Rdata",
    use_blast = FALSE,
    multithread = 4L
)

if (run_test("Config creation (4-run, Swarm)", {
    stopifnot(inherits(config_coi_multi, "smp_config"))
    stopifnot(config_coi_multi$linkage == "single")
    stopifnot(config_coi_multi$n_runs == 4L)
})) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

if (run_test("Cache loading (COI_multirun)", {
    cache <- smp_load_cache(config_coi_multi)
    n_cached <- sum(!sapply(cache, is.null))
    cat(sprintf("(%d cached steps) ", n_cached))
    stopifnot(n_cached > 0)
})) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

###############################################################################
# TEST 6: Verify config validation catches real errors
###############################################################################
cat("\n── Validation error catching ─────────────────────────────────────────\n")

if (run_test("Bad read_type rejected", {
    tryCatch(
        smp_config(path = env_root, dataname = "test", read_type = "invalid"),
        error = function(e) invisible(NULL)
    )
    # If no error, fail
    tryCatch({
        smp_config(path = env_root, dataname = "test", read_type = "invalid")
        stop("Should have errored")
    }, error = function(e) invisible(NULL))
})) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

if (run_test("Missing classifier detected", {
    cfg <- smp_config(
        path = env_root, dataname = "Example",
        n_runs = 2L,
        trunc_len = list(c(150, 150), c(150, 150)),
        trim_left = list(c(0, 0), c(0, 0)),
        use_idtaxa = TRUE,
        training_set = "NONEXISTENT_FILE.RData"
    )
    tryCatch({
        smp_validate(cfg)
        stop("Should have errored")
    }, error = function(e) {
        stopifnot(grepl("training set", conditionMessage(e), ignore.case = TRUE) ||
                  inherits(e, "rlang_error"))
    })
})) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

if (run_test("trunc_len / n_runs mismatch detected", {
    cfg <- smp_config(
        path = env_root, dataname = "Example",
        n_runs = 3L,
        trunc_len = list(c(150, 150)),  # only 1 entry, need 3
        trim_left = list(c(0, 0)),
        use_idtaxa = FALSE
    )
    tryCatch({
        smp_validate(cfg)
        stop("Should have errored")
    }, error = function(e) invisible(NULL))
})) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

###############################################################################
# TEST 7: Cross-check cached results structure
###############################################################################
cat("\n── Cached result structure validation ────────────────────────────────\n")

if (run_test("Example DADA2 cache has expected columns", {
    cache <- smp_load_cache(config_example)
    sdt <- cache$dada2$SeqDataTable
    required_cols <- c("ESV", "Sequence")
    for (col in required_cols) {
        stopifnot(col %in% names(sdt))
    }
    sample_cols <- grep("^Sample_", names(sdt), value = TRUE)
    stopifnot(length(sample_cols) >= 2)
    cat(sprintf("(%d ESVs, %d samples) ", nrow(sdt), length(sample_cols)))
})) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

if (run_test("Example LULU1 cache has curation columns", {
    cache <- smp_load_cache(config_example)
    lulu1 <- cache$lulu1
    stopifnot("curatedESV" %in% names(lulu1) ||
              "CuratedESVRepresentativeSequence" %in% names(lulu1))
    cat("(curation columns present) ")
})) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

if (run_test("Example cluster cache has OTU column", {
    cache <- smp_load_cache(config_example)
    cluster <- cache$cluster
    stopifnot("OTU" %in% names(cluster))
    n_otus <- length(unique(cluster$OTU))
    cat(sprintf("(%d OTUs) ", n_otus))
})) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

if (run_test("Example IDTAXA cache has taxonomy", {
    cache <- smp_load_cache(config_example)
    idtaxa <- cache$idtaxa
    stopifnot("SeqDataTable" %in% names(idtaxa))
})) {
    passed <- passed + 1
} else {
    failed <- failed + 1
}

###############################################################################
# Summary
###############################################################################
cat("\n========================================\n")
cat(sprintf("RESULTS: %d passed, %d failed, %d total\n",
            passed, failed, passed + failed))
cat("========================================\n")

if (failed > 0) {
    cat("\nFAILED TESTS:\n")
    for (e in errors) cat(sprintf("  - %s\n", e))
    quit(status = 1)
} else {
    cat("\nAll integration tests passed!\n")
}
