## meBART: Bayesian Additive Regression Trees with Measurement Error
## Copyright (C) 2025 Kevin McCoy, Zachary Wooten, and Christine Peterson

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/GPL-2




#' Fit a Bayesian Additive Regression Tree (BART) Model
#'
#' This function fits a Bayesian Additive Regression Tree model for continuous data using a set of explanatory variables.
#' It provides options for handling test data, priors, and other model parameters, including options for variable selection and sparsity control.
#'
#' @param x.train A matrix or data frame containing the explanatory variables for the training data.
#' @param y.train A numeric vector containing the continuous outcome variable for the training data.
#' @param x.test A matrix or data frame of test data (optional, default is an empty matrix). The structure should match `x.train`.
#' @param sparse A logical value (default is FALSE). If TRUE, uses a Dirichlet prior for variable selection.
#' @param theta A numeric value (default is 0). Parameter for model specification (unspecified usage).
#' @param omega A numeric value (default is 1). Parameter for model specification (unspecified usage).
#' @param a A numeric value (default is 0.5). Beta parameter used in DART (Dropout Additive Regression Trees).
#' @param b A numeric value (default is 1). Beta parameter used in DART.
#' @param augment A logical value (default is FALSE). An additional parameter (unspecified usage).
#' @param rho A numeric value (default is NULL). Sparsity parameter used in DART.
#' @param xinfo A matrix of cutpoints (default is an empty matrix). Specifies cutpoints for the model if desired.
#' @param usequants A logical value (default is FALSE). If TRUE, uniform quantiles are used for cutpoints; otherwise, uniform cutpoints are used.
#' @param cont A logical value (default is FALSE). Specifies if all variables are continuous.
#' @param rm.const A logical value (default is TRUE). If TRUE, removes constant columns from `x.train`.
#' @param sigest A numeric value (default is NA). Prior estimate for sigma; if NA, it is estimated from the data.
#' @param sigdf A numeric value (default is 3). Degrees of freedom for the sigma prior.
#' @param sigquant A numeric value (default is 0.90). Quantile for the sigma prior.
#' @param k A numeric value (default is 2.0). Prior parameter for the mean of the output data.
#' @param power A numeric value (default is 2.0). Tree depth prior, beta in the CGM (Conditional Gaussian Model).
#' @param base A numeric value (default is 0.95). Tree depth prior, alpha in the CGM.
#' @param sigmaf A numeric value (default is NA). Parameter for the sigma prior, not scaled by the number of trees.
#' @param lambda A numeric value (default is NA). Scale parameter used in the sigma prior.
#' @param fmean A numeric value (default is the mean of `y.train`). The mean of the output variable `y.train`.
#' @param w A numeric vector (default is a vector of ones). Weights to multiply the standard deviation.
#' @param ntree An integer (default is 200L). The number of trees to use in the model.
#' @param numcut An integer (default is 100L). The number of cutpoints to consider for each variable.
#' @param ndpost An integer (default is 1000L). The number of posterior draws to return.
#' @param nskip An integer (default is 100L). The number of observations for burn-in (MCMC).
#' @param keepevery An integer (default is 1L). The thinning parameter for MCMC.
#' @param nkeeptrain An integer (default is equal to `ndpost`). The number of training data posterior draws to keep.
#' @param nkeeptest An integer (default is equal to `ndpost`). The number of test data posterior draws to keep.
#' @param nkeeptestmean An integer (default is equal to `ndpost`). The number of MCMC iterations to return for the test mean.
#' @param nkeeptreedraws An integer (default is equal to `ndpost`). The number of MCMC iterations to return for tree draws.
#' @param printevery An integer (default is 100L). The number of iterations before progress is printed.
#' @param transposed A logical value (default is FALSE). Used if called by `mc.wbart`.
#'
#' @return A list with the following elements:
#' \item{mu}{The mean of the outcome variable `y.train`, centered back in the predictions.}
#' \item{yhat.train.mean}{The mean of the posterior predictions for the training data.}
#' \item{yhat.train}{The posterior predictions for the training data, centered.}
#' \item{yhat.test.mean}{The mean of the posterior predictions for the test data.}
#' \item{yhat.test}{The posterior predictions for the test data, centered.}
#' \item{varcount}{A matrix of variable counts in the model.}
#' \item{varprob}{A matrix of variable probabilities in the model.}
#' \item{treedraws}{A list of posterior tree draws (if `nkeeptreedraws > 0`).}
#' \item{x_draws}
#' \item{acceptances}
#' \item{proc.time}{The time taken for the C++ code to run.}
#' \item{varcount.mean}{The mean of the variable counts.}
#' \item{varprob.mean}{The mean of the variable probabilities.}
#' \item{rm.const}{A vector specifying which columns were removed from `x.train`.}
#' @export
#' @importFrom stats lm
#'
#' @examples
#' # Example usage:
#' set.seed(42)
#' x.train <- matrix(rnorm(1000), ncol = 10)
#' y.train <- rnorm(100)
#' result <- mebart(x.train, y.train, ntree = 50)
#' print(result$yhat.train.mean)
#'
mebart <- function(x.train, # explanatory variables for training data, matrix or dataframe
                   y.train, # continuous outcome variable
                   x.test = matrix(0.0, 0, 0), # x test data, same structure as x.train
                   sparse = FALSE, # use Dirichlet prior for variable selection
                   theta = 0, # TODO ???
                   omega = 1, # TODO ???
                   a = 0.5, # Beta parameters, used in DART only
                   b = 1, # Beta parameters, used in DART only
                   augment = FALSE, # TODO ???
                   rho = NULL, # Sparsity parameter, used in DART only
                   xinfo = matrix(0.0, 0, 0), # cutpoints, if you want to specify them
                   usequants = FALSE, # if false makes cutpoints uniform, if true uniform quantiles are used
                   cont = FALSE, # whether or not to assume all variables are continuous
                   rm.const = TRUE, # remove constant columns from x.train
                   sigest = NA, # sigma prior, if NA it is estimated from the data
                   sigdf = 3, # degrees of freedom for sigma prior
                   sigquant = .90, # quantile for sigma prior
                   k = 2.0, # k used in mu prior
                   power = 2.0, # tree depth prior, beta in CGM
                   base = .95, # tree depth prior, alpha in CGM
                   sigmaf = NA, # sigma from mu prior, not yet scaled by number of trees
                   lambda = NA, # scale parameter used in sigma prior
                   fmean = mean(y.train), # Mean of output data
                   w = rep(1, length(y.train)), # vector of weights to multiply the standard deviation, defaults to one
                   ntree = 200L, # number of trees
                   numcut = 100L, # number of cutpoints to consider
                   ndpost = 1000L, # number of posterior draws returned
                   nskip = 100L, # number of observations for burn-in
                   keepevery = 1L, # thinning parameter
                   nkeeptrain = ndpost, # number of training data posterior draws to keep
                   nkeeptest = ndpost, # number of test data posterior draws to keep
                   nkeeptestmean = ndpost, # number of MCMC iters to be returned for test mean
                   nkeeptreedraws = ndpost, # number of MCMC iters to be returned for tree draws
                   printevery = 100L, # print progress every printevery iterations
                   transposed = FALSE, # used if called by mc.wbart
                   proposal_sigma = meas_error_sigma, # standard deviation of the proposal distribution
                   meas_error_sigma, # standard deviation of the measurement error
                   x_mu,
                   x_sigma
                   ) {
    #--------------------------------------------------
    # data
    n <- length(y.train)

    if (!transposed) { # if being called wbart, skipped if being called by mc.wbart
        temp <- bartModelMatrix(
            x.train,
            numcut,
            usequants = usequants, # if false makes cutpoints uniform, if true uniform quantiles are used
            cont = cont, # whether or not to assume all variables are continuous
            xinfo = xinfo, # cutpoints, if you want to specify them
            rm.const = rm.const # remove constant columns from x.train
        )
        x.train <- t(temp$X) # transpose data matrix, has been converted to matrix
        numcut <- temp$numcut # number of cutpoints, vector of each variable
        xinfo <- temp$xinfo # cutpoints, now specified
        if (length(x.test) > 0) {
            x.test <- bartModelMatrix(x.test)
            x.test <- t(x.test[, temp$rm.const])
        }
        rm.const <- temp$rm.const # removed variables or all if kept
        grp <- temp$grp # counts out which variables are unique / which ones are factors
        rm(temp) # clears temp from R's workspace
    } else {
        rm.const <- NULL
        grp <- NULL
    }

    if (n != ncol(x.train)) { # check if the number of observations in x.train is equal to the length of y.train, remember x.train is transposed
        stop("The length of y.train and the number of rows in x.train must be identical")
    }

    p <- nrow(x.train)
    np <- ncol(x.test)
    if (length(rho) == 0) {
        rho <- p
    }
    if (length(rm.const) == 0) {
        rm.const <- 1:p
    }
    if (length(grp) == 0) {
        grp <- 1:p
    }

    # I added these
    if (nrow(meas_error_sigma) != ncol(meas_error_sigma)) {
        stop("The meas_error_sigma matrix must be square")
    }
    if (p != nrow(meas_error_sigma)){
        stop("The number of rows/columns in meas_error_sigma must be equal to the number of columns in x.train")
    }
    if (any(eigen(meas_error_sigma)$values <= 0) | any(meas_error_sigma != t(meas_error_sigma))) {
        stop("The meas_error_sigma matrix must be symmetric and positive definite")
    }


    y.train <- y.train - fmean # centers output data
    #------------------------------------------------------------------------------------------------------------
    # set nkeeps for thinning ######### Basically makes sure these variables dont conflict with one another
    if ((nkeeptrain != 0) & ((ndpost %% nkeeptrain) != 0)) { # if
        nkeeptrain <- ndpost
        cat("*****nkeeptrain set to ndpost\n") # cat is just a simpler version of print
    }
    if ((nkeeptest != 0) & ((ndpost %% nkeeptest) != 0)) {
        nkeeptest <- ndpost
        cat("*****nkeeptest set to ndpost\n")
    }
    if ((nkeeptestmean != 0) & ((ndpost %% nkeeptestmean) != 0)) {
        nkeeptestmean <- ndpost
        cat("*****nkeeptestmean set to ndpost\n")
    }
    if ((nkeeptreedraws != 0) & ((ndpost %% nkeeptreedraws) != 0)) {
        nkeeptreedraws <- ndpost
        cat("*****nkeeptreedraws set to ndpost\n")
    }
    #----------------------------------------------------------------------------------------------
    # prior
    nu <- sigdf # nu / DOF in sigma prior
    if (is.na(lambda)) { # NA by default
        if (is.na(sigest)) { # NA by default
            if (p < n) {
                df <- data.frame(t(x.train), y.train)
                lmf <- lm(y.train ~ ., df) # trains linear model
                sigest <- summary(lmf)$sigma # gets sigma from linear model
            } else {
                sigest <- sd(y.train)
            }
        }
        qchi <- qchisq(1.0 - sigquant, nu) # gets quantile from chi squared distribution
        lambda <- (sigest * sigest * qchi) / nu # lambda parameter for sigma prior
    } else {
        sigest <- sqrt(lambda)
    }

    if (is.na(sigmaf)) { # NA by default
        tau <- (max(y.train) - min(y.train)) / (2 * k * sqrt(ntree))
    } else {
        tau <- sigmaf / sqrt(ntree)
    }
    #-------------------------------------------------------------------------------------------------------
    ptm <- proc.time() # how much real and CPU time (in seconds) the currently running R process has already taken.
    # call
    res <- .Call(
        "cmebart",
        n, # number of observations in training data
        p, # dimension of x
        np, # number of observations in test data
        x.train, # pxn training data x
        y.train, # pxn training data x
        x.test, # p*np test data x
        ntree, # number of trees
        numcut, # number of cutpoints
        ndpost * keepevery, # total number of draws
        nskip, # number of observations for burn-in
        power, # tree depth prior, beta in CGM
        base, # tree depth prior, alpha in CGM
        tau, # sigma from mu prior, now scaled by number of trees
        nu, # DOF in sigma prior
        lambda, # scale parameter used in sigma prior
        sigest, # sigma prior, estimated from the data
        w, # vector of weights to multiply the standard deviation, defaults to one
        sparse, # use Dirichlet prior for variable selection
        theta, # TODO ???
        omega, # TODO ???
        grp, # counts out which variables are unique / which ones are factors
        a, # Beta parameters, used in DART only
        b, # Beta parameters, used in DART only
        rho, # Sparsity parameter, used in DART only
        augment, # TODO ???
        nkeeptrain, # number of training data posterior draws to keep
        nkeeptest, # number of test data posterior draws to keep
        nkeeptestmean, # number of MCMC iters to be returned for test mean
        nkeeptreedraws, # number of MCMC iters to be returned for tree draws
        printevery, # print progress every printevery iterations
        xinfo, # cutpoints, now specified
        proposal_sigma, # standard deviation of the proposal distribution
        meas_error_sigma, # standard deviation of the measurement error
        x_mu,
        x_sigma
    )

    res$proc.time <- proc.time() - ptm # how much time C++ code has taken

    # Centering variable, add back in
    res$mu <- fmean
    res$yhat.train.mean <- res$yhat.train.mean + fmean
    res$yhat.train <- res$yhat.train + fmean
    res$yhat.test.mean <- res$yhat.test.mean + fmean
    res$yhat.test <- res$yhat.test + fmean



    if (nkeeptreedraws > 0) {
        names(res$treedraws$cutpoints) <- dimnames(x.train)[[1]]
    }
    dimnames(res$varcount)[[2]] <- as.list(dimnames(x.train)[[1]])
    dimnames(res$varprob)[[2]] <- as.list(dimnames(x.train)[[1]])


    res$varcount.mean <- apply(res$varcount, 2, mean)
    res$varprob.mean <- apply(res$varprob, 2, mean)
    res$rm.const <- rm.const
    attr(res, "class") <- "mebart"
    return(res)
}
