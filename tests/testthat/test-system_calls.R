test_that("smp_system raises error on failed command", {
    expect_error(
        smp_system("false", character(0), "test_failure"),
        "failed"
    )
})

test_that("smp_system succeeds on valid command", {
    expect_no_error(
        smp_system("true", character(0), "test_success")
    )
})

test_that("smp_shell raises error on failed command", {
    expect_error(
        smp_shell("false", "test_failure"),
        "failed"
    )
})

test_that("smp_shell succeeds on valid command", {
    expect_no_error(
        smp_shell("true", "test_success")
    )
})
