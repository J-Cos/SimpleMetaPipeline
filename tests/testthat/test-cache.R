test_that("smp_cache and smp_load_cache round-trip correctly", {
    tmp <- tempdir()
    cfg <- smp_config(path = tmp, dataname = "cache_test")
    dir.create(file.path(tmp, "IntermediateOutputs"), showWarnings = FALSE)

    test_data <- list(a = 1, b = "hello", c = data.frame(x = 1:3))
    smp_cache(test_data, cfg, "test_step")

    # Check file exists
    expected_path <- file.path(
        tmp, "IntermediateOutputs", "cache_test_test_step.RDS"
    )
    expect_true(file.exists(expected_path))

    # Load and compare
    loaded <- readRDS(expected_path)
    expect_equal(loaded, test_data)

    # Clean up
    file.remove(expected_path)
})

test_that("smp_load_cache returns NULL for missing steps", {
    tmp <- tempdir()
    cfg <- smp_config(path = tmp, dataname = "empty_cache_test")
    dir.create(file.path(tmp, "IntermediateOutputs"), showWarnings = FALSE)

    cache <- smp_load_cache(cfg)
    expect_type(cache, "list")
    expect_true(all(sapply(cache, is.null)))
})

test_that("smp_clear_cache removes cached files", {
    tmp <- tempdir()
    cfg <- smp_config(path = tmp, dataname = "clear_test")
    dir.create(file.path(tmp, "IntermediateOutputs"), showWarnings = FALSE)

    # Create some fake caches
    smp_cache(list(x = 1), cfg, "dada2")
    smp_cache(list(y = 2), cfg, "lulu1")

    # Verify they exist
    expect_true(file.exists(
        file.path(tmp, "IntermediateOutputs", "clear_test_dada2.RDS")
    ))

    # Clear
    smp_clear_cache(cfg)

    # Verify gone
    expect_false(file.exists(
        file.path(tmp, "IntermediateOutputs", "clear_test_dada2.RDS")
    ))
    expect_false(file.exists(
        file.path(tmp, "IntermediateOutputs", "clear_test_lulu1.RDS")
    ))
})
