###########################
# SimpleMetaPipeline v3 — Control Script Template
###########################
#
# This template replaces the old ControlScript approach.
# All parameters are passed via smp_config() — no global variables.
#
# Usage:
#   conda activate ./BioinformaticPipeline/env
#   Rscript "BioinformaticPipeline/inst/templates/control_script_template.R"
#
###########################

library(SimpleMetaPipeline)

config <- smp_config(

    # ── General settings ──────────────────────────────────────────────────
    path     = "../BioinformaticPipeline_Env",
    dataname = "Example",
    # data_path = NULL,        # set if FASTQs are stored elsewhere
    multithread = 4L,

    # ── Cutadapt ──────────────────────────────────────────────────────────
    use_cutadapt = FALSE,
    fwd_primer   = "ACTG",
    rev_primer   = "ACTG",

    # ── DADA2 ─────────────────────────────────────────────────────────────
    n_runs    = 2L,
    trunc_len = list(c(150, 150), c(150, 150)),
    trim_left = list(c(0, 0), c(0, 0)),
    max_n     = 0L,
    max_ee    = c(2, 2),
    trunc_q   = 2L,
    desired_seq_length_range = NULL,
    pool      = "pseudo",
    mixed_orientation = FALSE,
    read_type = "paired-end",

    # ── LULU (stage 1 — ESV curation) ─────────────────────────────────────
    lulu1_match_rate           = 90,
    lulu1_min_rel_cooccurrence = 0.95,
    lulu1_ratio_type           = "min",

    # ── Clustering ────────────────────────────────────────────────────────
    linkage              = "complete",
    differences          = 1L,
    similarity_threshold = 0.97,
    threads              = 20L,

    # ── LULU (stage 2 — OTU curation) ─────────────────────────────────────
    lulu2_match_rate           = 90,
    lulu2_min_rel_cooccurrence = 0.95,
    lulu2_ratio_type           = "min",

    # ── IDTAXA ────────────────────────────────────────────────────────────
    use_idtaxa       = TRUE,
    training_set     = "GTDB_r207-mod_April2022.RData",
    desired_ranks    = c("rootrank", "domain", "phylum", "class", "order",
                         "family", "genus"),
    seqs_to_assign   = "ESVs",
    idtaxa_threshold = 40,
    idtaxa_parallel  = FALSE,

    # ── BLAST ─────────────────────────────────────────────────────────────
    use_blast        = FALSE,
    blast_db         = NULL,
    blast_thresholds = c(0, 0, 80, 85, 90, 95, 97),
    blast_desired_ranks = c("Root", "Kingdom", "Phylum", "Class", "Order",
                            "Family", "Genus", "Species")
)

# Run the pipeline
smp_run(config)
