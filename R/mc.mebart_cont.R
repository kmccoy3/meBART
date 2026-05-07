## meBART: Bayesian Additive Regression Trees with Measurement Error
## Copyright (C) 2025 Kevin McCoy, Zachary Wooten, and Christine Peterson

## This package is a modification of the BART package originally created by
## Robert McCulloch and Rodney Sparapani, who own the original copyright.
## Copyright (C) 2017 Robert McCulloch and Rodney Sparapani

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/GPL-2


#' @title Fit a meBART Model with \emph{Continuous} Output Data (Multi-core Version)
#'
#' @description This function fits a Bayesian Additive Regression Tree (BART) model which also also models measurement error in the explanatory variables. 
#' This function specifically handles \emph{continuous} output data.
#'
#' @param x.train A matrix or data frame containing the explanatory variables for the training data.
#' @param y.train An integer vector containing the binary outcome variable for the training data.
#' @param x.test A matrix or data frame of test data (optional, default is an empty matrix). The structure should match `x.train`.
#' @param sparse Whether to use Dirichlet additive regression trees (DART).
#' @param theta alpha parameter in Linero 2018 equation (2), governs DART prior (DART ONLY).
#' @param omega omega parameter for DART.
#' @param a Beta prior parameters, used in DART only.
#' @param b Beta prior parameters, used in DART only.
#' @param augment FALSE uses Assumption 2.2 in Linero 2018, TRUE uses assumption 2.1 (DART ONLY).
#' @param rho rho sparsity parameter, used in DART only.
#' @param xinfo A matrix of cutpoints (default is an empty matrix). Specifies cutpoints for the model if desired.
#' @param usequants If TRUE, uniform quantiles are used for cutpoints; otherwise, uniform cutpoints are used.
#' @param cont Specifies if all variables are continuous.
#' @param rm.const If TRUE, removes constant columns from `x.train`.
#' @param sigest Prior estimate for sigma; if NA, it is estimated from the data.
#' @param sigdf  Degrees of freedom for the sigma prior.
#' @param sigquant Quantile for the sigma prior.
#' @param k Prior parameter for the mean of the output data.
#' @param power Power parameter for the BART tree-depth prior (beta in original BART paper).
#' @param base Base parameter for the BART tree-depth prior (alpha in original BART paper).
#' @param sigmaf Parameter for the sigma prior, not scaled by the number of trees.
#' @param lambda Scale parameter used in the sigma prior.
#' @param fmean A numeric value (default is the mean of `y.train`). The mean of the output variable `y.train`.
#' @param w A numeric vector (default is a vector of ones). Weights to multiply the standard deviation. Must be of length n.
#' @param ntree The number of trees to use in the model.
#' @param numcut The number of cutpoints to consider for each variable.
#' @param ndpost The number of posterior draws to return.
#' @param nskip The number of observations for burn-in (MCMC).
#' @param keepevery The thinning parameter for MCMC.
#' @param printevery How often to print progress during MCMC.
#' @param keeptrainfits Whether to keep the training fits in the output.
#' @param transposed Is the input data transposed? Used if called by `mc.wbart`.
#' @param mc.cores Number of cores to use for parallel processing.
#' @param nice Nice level for the process.
#' @param seed Random seed for reproducibility. 
#' @param proposal_sigma A numeric value (default is `meas_error_sigma`). Covariance matrix of the proposal distribution.
#' @param meas_error_sigma A numeric matrix. Covariance matrix of the measurement error, must be square and match the number of columns in `x.train`.
#' @param x_mu A numeric column vector. Mean of the explanatory variables.
#' @param x_sigma A numeric matrix. Covariance matrix of the explanatory variables.
#'
#' @return A list with the following elements:
#' \item{mu}{The mean of the outcome variable `y.train`.}
#' \item{yhat.train.mean}{The mean of the posterior predictions for the training data.}
#' \item{yhat.train}{The posterior predictions for the training data.}
#' \item{yhat.test.mean}{The mean of the posterior predictions for the test data.}
#' \item{yhat.test}{The posterior predictions for the test data.}
#' \item{varcount}{An (ndpost * p) matrix of variable counts in the model.}
#' \item{varprob}{An (ndpost * p) matrix of variable probabilities in the model.}
#' \item{treedraws}{A list of posterior tree draws.}
#' \item{x_draws}{Tensor with posterior draws for the training data.}
#' \item{acceptances}{Whether or not the proposal was accepted.}
#' \item{proc.time}{The time taken for the C++ code to run.}
#' \item{varcount.mean}{The mean of the variable counts.}
#' \item{varprob.mean}{The mean of the variable probabilities.}
#' \item{rm.const}{A vector specifying which columns were removed from `x.train`.}
#' \item{ndpost}{The number of posterior draws returned.}
#' @export
#' 
#' @importFrom parallel detectCores
#' @importFrom tools psnice
#' @importFrom abind abind
#'
#' @examples
#' # Example usage:
#' set.seed(0)
#' x.train <- matrix(rnorm(1000), ncol = 10)
#' y.train <- rnorm(100)
#' mdl <- mc.mebart_cont(x.train, y.train,
#'               ndpost = 500,
#'               ntree = 50,
#'               meas_error_sigma = diag(10),
#'               x_mu = matrix(0, nrow=10, ncol=1),
#'               x_sigma = diag(10))
#' preds <- mdl$yhat.train.mean
#' print(paste0("Training MSE: ", round(mean((y.train - preds)^2), 3)))
mc.mebart_cont <- function(
    x.train, 
    y.train, 
    x.test = matrix(0.0, 0, 0),
    sparse = FALSE, 
    theta = 0, 
    omega = 1,
    a = 0.5, 
    b = 1, 
    augment = FALSE, 
    rho = NULL,
    xinfo = matrix(0.0, 0, 0), 
    usequants = FALSE,
    cont = FALSE, 
    rm.const = TRUE,
    sigest = NA, 
    sigdf = 3, 
    sigquant = 0.90,
    k = 2.0, 
    power = 2.0, 
    base = .95,
    sigmaf = NA, 
    lambda = NA,
    fmean = mean(y.train),
    w = rep(1, length(y.train)),
    ntree = 200L, 
    numcut = 100L,
    ndpost = 1000L, 
    nskip = 100L,
    keepevery = 1L, 
    printevery = 100L,
    keeptrainfits = TRUE, 
    transposed = FALSE,
    mc.cores = 2L, 
    nice = 19L,
    seed = 99L,
    proposal_sigma = meas_error_sigma, # standard deviation of the proposal distribution
    meas_error_sigma, # standard deviation of the measurement error
    x_mu,
    x_sigma
    ) {
    if (.Platform$OS.type != "unix") {
        stop("parallel::mcparallel/mccollect do not exist on windows")
    }

    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()

    if (!transposed) {
        temp <- bartModelMatrix(x.train, numcut,
            usequants = usequants,
            cont = cont, xinfo = xinfo, rm.const = rm.const
        )
        x.train <- t(temp$X)
        numcut <- temp$numcut
        xinfo <- temp$xinfo
        ## if(length(x.test)>0)
        ##     x.test = t(bartModelMatrix(x.test[ , temp$rm.const]))
        if (length(x.test) > 0) {
            x.test <- bartModelMatrix(x.test)
            x.test <- t(x.test[, temp$rm.const])
        }
        rm.const <- temp$rm.const
        rm(temp)
    }

    mc.cores.detected <- detectCores()

    if (mc.cores > mc.cores.detected) mc.cores <- mc.cores.detected

    mc.ndpost <- ceiling(ndpost / mc.cores)

    for (i in 1:mc.cores) {
        parallel::mcparallel(
            {
                psnice(value = nice)
                mebart_cont(
                    x.train = x.train, y.train = y.train, x.test = x.test,
                    sparse = sparse, theta = theta, omega = omega,
                    a = a, b = b, augment = augment, rho = rho,
                    xinfo = xinfo,
                    sigest = sigest, sigdf = sigdf, sigquant = sigquant,
                    k = k, power = power, base = base,
                    sigmaf = sigmaf, lambda = lambda, fmean = fmean, w = w,
                    ntree = ntree, numcut = numcut,
                    ndpost = mc.ndpost, nskip = nskip, keepevery = keepevery,
                    printevery = printevery, transposed = TRUE, 
                    proposal_sigma = meas_error_sigma, # standard deviation of the proposal distribution
                   meas_error_sigma = meas_error_sigma, # standard deviation of the measurement error
                   x_mu = x_mu,
                   x_sigma = x_sigma
                )
            },
            ## treesaslists=treesaslists)},
            silent = (i != 1)
        )
        ## to avoid duplication of output
        ## capture stdout from first posterior only
    }

    post.list <- parallel::mccollect()

    post <- post.list[[1]]



    if (mc.cores == 1 | attr(post, "class") != "mebart_cont") {
        return(post)
    } else {
        if (class(rm.const)[1] != "logical") post$rm.const <- rm.const
        post$ndpost <- mc.cores * mc.ndpost
        p <- nrow(x.train[post$rm.const, ])


        old.text <- paste0(
            as.character(mc.ndpost), " ", as.character(ntree),
            " ", as.character(p)
        )
        old.stop <- nchar(old.text)

        post$treedraws$trees <- sub(
            old.text,
            paste0(
                as.character(post$ndpost), " ",
                as.character(ntree), " ",
                as.character(p)
            ),
            post$treedraws$trees
        )

        for (i in 2:mc.cores) {
            post$yhat.train <- rbind(
                post$yhat.train,
                post.list[[i]]$yhat.train
            )

            if (length(post$yhat.test) > 0) {
                post$yhat.test <- rbind(
                    post$yhat.test,
                    post.list[[i]]$yhat.test
                )
            }

            post$x_draws <- abind(
                post$x_draws,
                post.list[[i]]$x_draws,
                along=3
            )


            post$sigma <- cbind(post$sigma, post.list[[i]]$sigma)

            post$treedraws$trees <- paste0(
                post$treedraws$trees,
                substr(
                    post.list[[i]]$treedraws$trees, old.stop + 2,
                    nchar(post.list[[i]]$treedraws$trees)
                )
            )

            if (length(post$varcount) > 0) {
                post$varcount <- rbind(post$varcount, post.list[[i]]$varcount)
                post$varprob <- rbind(post$varprob, post.list[[i]]$varprob)
            }

            post$proc.time["elapsed"] <- max(
                post$proc.time["elapsed"],
                post.list[[i]]$proc.time["elapsed"]
            )
            for (j in 1:5) {
                if (j != 3) {
                    post$proc.time[j] <- post$proc.time[j] + post.list[[i]]$proc.time[j]
                }
            }
        }

        if (length(post$yhat.train.mean) > 0) {
            post$yhat.train.mean <- apply(post$yhat.train, 2, mean)
        }

        if (length(post$yhat.test.mean) > 0) {
            post$yhat.test.mean <- apply(post$yhat.test, 2, mean)
        }

        if (length(post$varcount) > 0) {
            post$varcount.mean <- apply(post$varcount, 2, mean)
            post$varprob.mean <- apply(post$varprob, 2, mean)
        }

        attr(post, "class") <- "mebart_cont"

        return(post)
    }
}
