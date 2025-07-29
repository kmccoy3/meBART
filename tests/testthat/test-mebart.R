

############################## Initial Test ##############################


test_that("test_that framework is working", {
    expect_identical(2 + 3, 5)
})


############################## 1-D meBART Test ##############################


test_that("1-D meBART objects are identical to previous commits", {
    
    # Pass output of functions to sink so they don't clutter the Test console
    f = file()
    sink(file = f)
    
    # Generate data
    set.seed(0)
    n <- 20
    data <- data.frame(X1 = rnorm(n), y = rnorm(n))
    
    # Run meBART model
    Sigma <- matrix(c(1), nrow = 1)
    set.seed(0)
    meBART_mdl_1D <- meBART::mebart(data[, 1], data[, 2], 
                                    meas_error_sigma = Sigma,
                                    x_mu = matrix(0),
                                    x_sigma = matrix(1))
    
    # Close sink
    closeAllConnections()
    
    # Load old model
    old_mdl_1D <- readRDS("meBART_mdl_1D.rds")
    
    # Run tests
    expect_equal(meBART_mdl_1D$sigma, old_mdl_1D$sigma)
    expect_equal(meBART_mdl_1D$yhat.train.mean, old_mdl_1D$yhat.train.mean)
    expect_equal(meBART_mdl_1D$yhat.train, old_mdl_1D$yhat.train)
    expect_equal(meBART_mdl_1D$varcount, old_mdl_1D$varcount)
    expect_equal(meBART_mdl_1D$varprob, old_mdl_1D$varprob)
    # expect_equal(meBART_mdl_1D$treedraws, old_mdl_1D$treedraws)
    expect_equal(meBART_mdl_1D$x_draws, old_mdl_1D$x_draws)
    expect_equal(meBART_mdl_1D$acceptances, old_mdl_1D$acceptances)
    expect_equal(meBART_mdl_1D$mu, old_mdl_1D$mu)
    expect_equal(meBART_mdl_1D$varcount.mean, old_mdl_1D$varcount.mean)
    expect_equal(meBART_mdl_1D$varprob.mean, old_mdl_1D$varprob.mean)
    expect_equal(meBART_mdl_1D$rm.const, old_mdl_1D$rm.const)
})


############################## 2-D meBART Test ##############################


test_that("2-D meBART objects are identical to previous commits", {
    
    # Pass output of functions to sink so they don't clutter the Test console
    f = file()
    sink(file = f)
    
    # Generate data
    set.seed(0)
    n <- 20
    data <- data.frame(X1 = rnorm(n), X2 = rnorm(n), y = rnorm(n))
    
    # Run meBART model
    Sigma <- matrix(c(1, 0, 0, 1), nrow = 2)
    set.seed(0)
    meBART_mdl_2D <- meBART::mebart(data[, 1:2], data[, 3], 
                                    meas_error_sigma = Sigma,
                                    x_mu = matrix(c(0, 0), nrow=2),
                                    x_sigma = matrix(c(1, 0, 0, 1), nrow=2))
    
    # Close sink
    closeAllConnections()
    
    # Load old model
    old_mdl_2D <- readRDS("meBART_mdl_2D.rds")
    
    # Run tests
    expect_equal(meBART_mdl_2D$sigma, old_mdl_2D$sigma)
    expect_equal(meBART_mdl_2D$yhat.train.mean, old_mdl_2D$yhat.train.mean)
    expect_equal(meBART_mdl_2D$yhat.train, old_mdl_2D$yhat.train)
    expect_equal(meBART_mdl_2D$varcount, old_mdl_2D$varcount)
    expect_equal(meBART_mdl_2D$varprob, old_mdl_2D$varprob)
    # expect_equal(meBART_mdl_2D$treedraws, old_mdl_2D$treedraws)
    expect_equal(meBART_mdl_2D$x_draws, old_mdl_2D$x_draws)
    expect_equal(meBART_mdl_2D$acceptances, old_mdl_2D$acceptances)
    expect_equal(meBART_mdl_2D$mu, old_mdl_2D$mu)
    expect_equal(meBART_mdl_2D$varcount.mean, old_mdl_2D$varcount.mean)
    expect_equal(meBART_mdl_2D$varprob.mean, old_mdl_2D$varprob.mean)
    expect_equal(meBART_mdl_2D$rm.const, old_mdl_2D$rm.const)
})

