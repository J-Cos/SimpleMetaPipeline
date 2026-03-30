test_that("resolve_cores returns detectCores for TRUE", {
    result <- resolve_cores(TRUE)
    expect_type(result, "integer")
    expect_true(result >= 1)
})

test_that("resolve_cores returns integer for numeric input", {
    expect_equal(resolve_cores(8), 8L)
    expect_equal(resolve_cores(1), 1L)
})

test_that("allOrients returns 4 orientations", {
    skip_if_not_installed("Biostrings")
    orients <- allOrients("ACGT")
    expect_length(orients, 4)
    expect_named(orients, c("Forward", "Complement", "Reverse", "RevComp"))
    expect_equal(orients[["Forward"]], "ACGT")
    expect_equal(orients[["RevComp"]], "ACGT")  # palindrome
})

test_that("allOrients handles non-palindromic primer", {
    skip_if_not_installed("Biostrings")
    orients <- allOrients("AACG")
    expect_equal(orients[["Forward"]], "AACG")
    expect_equal(orients[["RevComp"]], "CGTT")
})
