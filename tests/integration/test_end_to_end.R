#!/usr/bin/env Rscript
###############################################################################
# End-to-end test: Run the new pipeline on 18s_test from raw FASTQs,
# then compare with old cached results.
###############################################################################

devtools::load_all(".")

env_root <- normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
cat("Environment root:", env_root, "\n")

# ── Step 1: Load old results for comparison ──────────────────────────────────
cat("\n=== Loading old cached results ===\n")
old_dada <- readRDS(file.path(env_root, "IntermediateOutputs",
                               "18s_test_DadaOutput.RDS"))
old_lulu1 <- readRDS(file.path(env_root, "IntermediateOutputs",
                                "18s_test_LuluOutput1.RDS"))
old_cluster <- readRDS(file.path(env_root, "IntermediateOutputs",
                                  "18s_test_ClusterOutput.RDS"))
old_lulu2 <- readRDS(file.path(env_root, "IntermediateOutputs",
                                "18s_test_LuluOutput2.RDS"))
old_idtaxa <- readRDS(file.path(env_root, "IntermediateOutputs",
                                 "18s_test_IdtaxaOutput.RDS"))

cat("Old DadaOutput: ", class(old_dada), "\n")
if (is.data.frame(old_dada)) {
    cat("  dim:", nrow(old_dada), "x", ncol(old_dada), "\n")
    cat("  cols:", paste(head(names(old_dada), 6), collapse = ", "), "\n")
} else {
    cat("  names:", paste(names(old_dada), collapse = ", "), "\n")
}
cat("Old LuluOutput1:", nrow(old_lulu1), "x", ncol(old_lulu1), "\n")
cat("  cols:", paste(head(names(old_lulu1), 8), collapse = ", "), "\n")
cat("Old ClusterOutput:", nrow(old_cluster), "x", ncol(old_cluster), "\n")
cat("  cols:", paste(head(names(old_cluster), 8), collapse = ", "), "\n")
cat("Old LuluOutput2:", nrow(old_lulu2), "x", ncol(old_lulu2), "\n")
cat("Old IdtaxaOutput:", class(old_idtaxa), "\n")
if (is.data.frame(old_idtaxa)) {
    cat("  dim:", nrow(old_idtaxa), "x", ncol(old_idtaxa), "\n")
} else {
    cat("  names:", paste(names(old_idtaxa), collapse = ", "), "\n")
}

# ── Step 2: Create config and clear new cache ────────────────────────────────
cat("\n=== Creating config ===\n")
config <- smp_config(
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
    desired_ranks = c("rootrank", "kingdom", "division", "phylum", "class",
                      "order", "family", "genus"),
    seqs_to_assign = "ESVs",
    idtaxa_threshold = 40,
    use_blast = FALSE,
    multithread = 4L,
    seed = 42L
)
print(config)

# ── Step 3: Validate ─────────────────────────────────────────────────────────
cat("\n=== Validating ===\n")
smp_validate(config)

# ── Step 4: Clear new-format cache (keep legacy files intact) ────────────────
cat("\n=== Clearing new-format cache ===\n")
new_cache_files <- list.files(
    file.path(env_root, "IntermediateOutputs"),
    pattern = "^18s_test_(cutadapt|dada2|lulu1|cluster|lulu2|idtaxa|blast)\\.RDS$",
    full.names = TRUE
)
if (length(new_cache_files) > 0) {
    file.remove(new_cache_files)
    cat("Removed", length(new_cache_files), "new-format cache files\n")
}

# ── Step 5: Run the pipeline ─────────────────────────────────────────────────
cat("\n=== Running pipeline ===\n")
result <- smp_run(config)

# ── Step 6: Compare results ──────────────────────────────────────────────────
cat("\n=== Comparing results ===\n")

# Extract new DADA2 SeqDataTable — handle both new (list) and legacy (df)
if (is.data.frame(result$dada2)) {
    new_dada <- result$dada2
} else {
    new_dada <- result$dada2$SeqDataTable
}

# old_dada is from the legacy cache — it's a data frame directly
if (!is.data.frame(old_dada)) {
    if ("SeqDataTable" %in% names(old_dada)) {
        old_dada <- old_dada$SeqDataTable
    }
}

cat("New DADA2 result:", nrow(new_dada), "ESVs,", ncol(new_dada), "cols\n")
cat("Old DADA2 result:", nrow(old_dada), "ESVs,", ncol(old_dada), "cols\n")

if (nrow(new_dada) == nrow(old_dada)) {
    cat("  ✓ ESV count matches!\n")
} else {
    cat("  ✗ ESV count differs\n")
}

# Check sequences match
old_seqs <- sort(as.character(old_dada$Sequence))
new_seqs <- sort(as.character(new_dada$Sequence))
seq_match <- identical(old_seqs, new_seqs)
cat("  Sequences identical:", seq_match, "\n")

# Check sample columns
old_sample_cols <- sort(grep("^Sample_", names(old_dada), value = TRUE))
new_sample_cols <- sort(grep("^Sample_", names(new_dada), value = TRUE))
cat("  Old sample columns:", length(old_sample_cols), "\n")
cat("  New sample columns:", length(new_sample_cols), "\n")
cat("  Sample columns identical:", identical(old_sample_cols, new_sample_cols), "\n")

# Check LULU1
if (is.data.frame(result$lulu1)) {
    new_lulu1 <- result$lulu1
} else if ("SeqDataTable" %in% names(result$lulu1)) {
    new_lulu1 <- result$lulu1$SeqDataTable
} else {
    new_lulu1 <- result$lulu1
}
cat("\nNew LULU1 result:", nrow(new_lulu1), "rows,", ncol(new_lulu1), "cols\n")
cat("Old LULU1 result:", nrow(old_lulu1), "rows,", ncol(old_lulu1), "cols\n")
if (nrow(new_lulu1) == nrow(old_lulu1)) cat("  ✓ Row count matches!\n")

# Check clustering
if (is.data.frame(result$cluster)) {
    new_cluster <- result$cluster
} else {
    new_cluster <- result$cluster
}
cat("\nNew Cluster result:", nrow(new_cluster), "rows,",
    length(unique(new_cluster$OTU)), "OTUs\n")
cat("Old Cluster result:", nrow(old_cluster), "rows,",
    length(unique(old_cluster$OTU)), "OTUs\n")
if (length(unique(new_cluster$OTU)) == length(unique(old_cluster$OTU))) {
    cat("  ✓ OTU count matches!\n")
}

# Check LULU2
if (is.data.frame(result$lulu2)) {
    new_lulu2 <- result$lulu2
} else {
    new_lulu2 <- result$lulu2
}
cat("\nNew LULU2 result:", nrow(new_lulu2), "rows,", ncol(new_lulu2), "cols\n")
cat("Old LULU2 result:", nrow(old_lulu2), "rows,", ncol(old_lulu2), "cols\n")
if (nrow(new_lulu2) == nrow(old_lulu2)) cat("  ✓ Row count matches!\n")

# Check final output
cat("\n=== Final output ===\n")
new_final <- readRDS(file.path(env_root, "Results", "18s_test_SeqDataTable.RDS"))
cat("Final SeqDataTable:", nrow(new_final), "rows,", ncol(new_final), "cols\n")

cat("\n=== End-to-end test complete ===\n")
