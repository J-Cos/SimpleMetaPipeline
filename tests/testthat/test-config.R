test_that("smp_config creates valid config with defaults", {
    cfg <- smp_config(
        path = tempdir(),
        dataname = "test"
    )
    expect_s3_class(cfg, "smp_config")
    expect_equal(cfg$dataname, "test")
    expect_equal(cfg$multithread, 4L)
    expect_equal(cfg$seed, 42L)
    expect_equal(cfg$read_type, "paired-end")
    expect_equal(cfg$linkage, "complete")
    expect_false(cfg$use_cutadapt)
    expect_true(cfg$use_idtaxa)
    expect_false(cfg$use_blast)
    expect_equal(cfg$lulu1_ratio_type, "min")
    expect_equal(cfg$lulu2_ratio_type, "min")
    expect_equal(cfg$pool, "pseudo")
})

test_that("smp_config match.arg catches invalid read_type", {
    expect_error(
        smp_config(path = tempdir(), dataname = "test", read_type = "invalid"),
        "arg"
    )
})

test_that("smp_config match.arg catches invalid linkage", {
    expect_error(
        smp_config(path = tempdir(), dataname = "test", linkage = "average"),
        "arg"
    )
})

test_that("smp_config match.arg catches invalid seqs_to_assign", {
    expect_error(
        smp_config(path = tempdir(), dataname = "test",
                   seqs_to_assign = "contigs"),
        "arg"
    )
})

test_that("smp_config resolves data_path to path by default", {
    cfg <- smp_config(path = tempdir(), dataname = "test")
    expect_equal(cfg$data_path, normalizePath(tempdir(), mustWork = FALSE))
})

test_that("smp_config accepts explicit data_path", {
    dp <- file.path(tempdir(), "custom")
    cfg <- smp_config(path = tempdir(), dataname = "test", data_path = dp)
    expect_equal(cfg$data_path, normalizePath(dp, mustWork = FALSE))
})

test_that("smp_config coerces multithread TRUE to integer", {
    cfg <- smp_config(path = tempdir(), dataname = "test", multithread = TRUE)
    expect_type(cfg$multithread, "integer")
    expect_true(cfg$multithread >= 1)
})

test_that("smp_config coerces integer parameters", {
    cfg <- smp_config(path = tempdir(), dataname = "test",
                      seed = 99, n_runs = 2, differences = 3)
    expect_type(cfg$seed, "integer")
    expect_equal(cfg$seed, 99L)
    expect_type(cfg$n_runs, "integer")
    expect_equal(cfg$n_runs, 2L)
    expect_type(cfg$differences, "integer")
    expect_equal(cfg$differences, 3L)
})

test_that("print.smp_config runs without error", {
    cfg <- smp_config(path = tempdir(), dataname = "test")
    expect_no_error(capture.output(print(cfg), type = "message"))
})
