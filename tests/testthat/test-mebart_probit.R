


############################## pBART Test ##############################

test_that("pbart objects are identical to previous commits", {
    
    # Pass output of functions to sink so they don't clutter the Test console
    f = file()
    sink(file = f)
    
    # Generate data
    set.seed(0)
    n <- 20
    data <- data.frame(X1 = rnorm(n), X2 = rnorm(n), y = rbinom(n, 1, 0.5))
    
    # Run meBART model
    Sigma <- matrix(c(1, 0, 0, 1), nrow = 2)
    set.seed(0)
    pbart_mdl <- meBART::mebart_probit(data[, 1:2], data[, 3], 
                                    meas_error_sigma = Sigma,
                                    x_mu = matrix(c(0, 0), nrow=2),
                                    x_sigma = matrix(c(1, 0, 0, 1), nrow=2))
    
    # Close sink
    closeAllConnections()
    
    # Load old model
    old_mdl <- readRDS("pbart_mdl.rds")
    
    # Run tests
    expect_equal(pbart_mdl$yhat.train, old_mdl$yhat.train)
    expect_equal(pbart_mdl$varcount, old_mdl$varcount)
    expect_equal(pbart_mdl$varprob, old_mdl$varprob)
    # expect_equal(pbart_mdl$treedraws, old_mdl$treedraws)
    expect_equal(pbart_mdl$x_draws, old_mdl$x_draws)
    expect_equal(pbart_mdl$acceptances, old_mdl$acceptances)
    expect_equal(pbart_mdl$prob.train, old_mdl$prob.train)
    expect_equal(pbart_mdl$prob.train.mean, old_mdl$prob.train.mean)
    expect_equal(pbart_mdl$varcount.mean, old_mdl$varcount.mean)
    expect_equal(pbart_mdl$varprob.mean, old_mdl$varprob.mean)
    expect_equal(pbart_mdl$rm.const, old_mdl$rm.const)
    expect_equal(pbart_mdl$binaryOffset, old_mdl$binaryOffset)
})

