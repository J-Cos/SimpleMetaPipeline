test_that("smp_validate catches missing FASTQ directory", {
    cfg <- smp_config(
        path = tempdir(),
        dataname = "nonexistent_dataset"
    )
    dir.create(file.path(tempdir(), "IntermediateOutputs"), showWarnings = FALSE)
    dir.create(file.path(tempdir(), "Results"), showWarnings = FALSE)

    expect_error(smp_validate(cfg), class = "rlang_error")
})

test_that("smp_validate catches trunc_len / n_runs mismatch", {
    cfg <- smp_config(
        path = tempdir(),
        dataname = "test",
        n_runs = 3L,
        trunc_len = list(c(150, 150))  # only 1, not 3
    )
    dir.create(file.path(tempdir(), "IntermediateOutputs"), showWarnings = FALSE)
    dir.create(file.path(tempdir(), "Results"), showWarnings = FALSE)

    expect_error(smp_validate(cfg), class = "rlang_error")
})

test_that("smp_validate catches missing IDTAXA training set", {
    tmp <- tempdir()
    cfg <- smp_config(
        path = tmp,
        dataname = "test",
        use_idtaxa = TRUE,
        training_set = "nonexistent.RData"
    )
    dir.create(file.path(tmp, "FASTQs", "test", "Run1"), recursive = TRUE,
               showWarnings = FALSE)
    dir.create(file.path(tmp, "IntermediateOutputs"), showWarnings = FALSE)
    dir.create(file.path(tmp, "Results"), showWarnings = FALSE)
    dir.create(file.path(tmp, "Data", "Classifiers"), recursive = TRUE,
               showWarnings = FALSE)

    expect_error(smp_validate(cfg), class = "rlang_error")
})

test_that("smp_validate catches missing BLAST database", {
    tmp <- tempdir()
    cfg <- smp_config(
        path = tmp,
        dataname = "test",
        use_idtaxa = FALSE,
        use_blast = TRUE,
        blast_db = "nonexistent_db"
    )
    dir.create(file.path(tmp, "FASTQs", "test", "Run1"), recursive = TRUE,
               showWarnings = FALSE)
    dir.create(file.path(tmp, "IntermediateOutputs"), showWarnings = FALSE)
    dir.create(file.path(tmp, "Results"), showWarnings = FALSE)

    expect_error(smp_validate(cfg), class = "rlang_error")
})

test_that("smp_validate passes on valid configuration", {
    tmp <- tempdir()
    cfg <- smp_config(
        path = tmp,
        dataname = "valid_test",
        use_idtaxa = FALSE,
        use_blast = FALSE,
        n_runs = 1L,
        trunc_len = list(c(150, 150)),
        trim_left = list(c(0, 0))
    )
    dir.create(file.path(tmp, "FASTQs", "valid_test", "Run1"),
               recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(tmp, "IntermediateOutputs"), showWarnings = FALSE)
    dir.create(file.path(tmp, "Results"), showWarnings = FALSE)

    expect_no_error(smp_validate(cfg))
})
