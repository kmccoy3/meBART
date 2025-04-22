
########## Initial Test ##########

test_that("test_that framework is working", {
    expect_identical(2+3, 5)
})

########## meBART Test ##########

# f = file()
# sink(file=f)
# 
# # Fix random seed to make the results reproducible
# set.seed(0)
# 
# # number of observations
# n <- 50
# 
# # variance of measurement error and output, respectively
# sigma_meas <- 0.1
# sigma_y <- 0.1
# 
# # Generate x values uniformly and add measurement error
# x_true <- runif(n, 0, 1)
# x_meas <- x_true + rnorm(n, 0, sigma_meas)
# y_true <- 1 * (x_true >= 0.5)
# y_obs <- y_true + rnorm(n, 0, sigma_y)
# 
# # Convert to df
# df_me <- data.frame(X=x_meas, y=y_obs)
# # df_no_me <- data.frame(X=x_true, y=y_obs)
# 
# # Parameters
# ndpost <- 500 # (Default: 1000L)
# nskip <- 100 # (Default: 100L)
# ntree <- 100 # (Default: 200L)
# 
# meas_error_sd = 0.10
# 
# # Data WITH Measurement Error
# set.seed(0)
# me_mdl <- meBART::mebart(df_me$X, df_me$y, 
#                          nskip=nskip, ndpost=ndpost, ntree=ntree, 
#                          meas_error_sd=meas_error_sd)
# 
# old_mdl <- readRDS("me_mdl.rds")
# 
# sink()
# 
# 
# test_that("meBART objects are identical to previous commits", {
#     expect_equal(me_mdl$sigma, old_mdl$sigma)
#     expect_equal(me_mdl$yhat.train.mean, old_mdl$yhat.train.mean)
#     expect_equal(me_mdl$yhat.train, old_mdl$yhat.train)
#     expect_equal(me_mdl$varcount, old_mdl$varcount)
#     expect_equal(me_mdl$varprob, old_mdl$varprob)
#     # expect_equal(me_mdl$treedraws, old_mdl$treedraws)
#     expect_equal(me_mdl$x_draws, old_mdl$x_draws)
#     expect_equal(me_mdl$acceptances, old_mdl$acceptances)
#     expect_equal(me_mdl$mu, old_mdl$mu)
#     expect_equal(me_mdl$varcount.mean, old_mdl$varcount.mean)
#     expect_equal(me_mdl$varprob.mean, old_mdl$varprob.mean)
#     expect_equal(me_mdl$rm.const, old_mdl$rm.const)
# })

########## Multivariate meBART Test ##########

f = file()
sink(file=f)

set.seed(0)
n <- 20
data <- data.frame(X1=rnorm(n), X2=rnorm(n), y= rnorm(n))

set.seed(0)
Sigma <- matrix(c(1, 0, 0, 1), nrow=2)
me_mdl <- meBART::mebart(data[,1:2], data[,3], meas_error_sigma=Sigma)

closeAllConnections()

old_mdl <- readRDS("me_multivariate_mdl.rds")



test_that("meBART objects are identical to previous commits", {
    expect_equal(me_mdl$sigma, old_mdl$sigma)
    expect_equal(me_mdl$yhat.train.mean, old_mdl$yhat.train.mean)
    expect_equal(me_mdl$yhat.train, old_mdl$yhat.train)
    expect_equal(me_mdl$varcount, old_mdl$varcount)
    expect_equal(me_mdl$varprob, old_mdl$varprob)
    # expect_equal(me_mdl$treedraws, old_mdl$treedraws)
    expect_equal(me_mdl$x_draws, old_mdl$x_draws)
    expect_equal(me_mdl$acceptances, old_mdl$acceptances)
    expect_equal(me_mdl$mu, old_mdl$mu)
    expect_equal(me_mdl$varcount.mean, old_mdl$varcount.mean)
    expect_equal(me_mdl$varprob.mean, old_mdl$varprob.mean)
    expect_equal(me_mdl$rm.const, old_mdl$rm.const)
})

