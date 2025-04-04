
# Test self-defined functions in mefuns.cpp file

###################################################################
############################## cmin ##############################

cmin_wrapper <- function(a, b){
    return(.Call("cmin", a, b, PACKAGE="meBART"))
}

test_that("min function is working", {
    expect_identical(cmin_wrapper(0, 0), 0)
    expect_identical(cmin_wrapper(1, 0), 0)
    expect_identical(cmin_wrapper(0, 1), 0)
    expect_identical(cmin_wrapper(-1, 0), -1)
    expect_identical(cmin_wrapper(-1, 1), -1)
    expect_identical(cmin_wrapper(1, -1), -1)
})

###################################################################
############################## Hello ##############################

dnorm_wrapper <- function(x, mean, sd){
    return(.Call("cdnorm", x, mean, sd, PACKAGE="meBART"))
}

test_that("dnorm function is working", {
    expect_equal(dnorm_wrapper(0, 0, 1), 0.39894228)
    # expect_equal(dnorm_wrapper(0, 0, -1), 0.39894228)
    # expect_error(dnorm_wrapper(1, 1, 0), "Sigma cannot be 0!")
})
