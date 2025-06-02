## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017 Robert McCulloch and Rodney Sparapani

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


#' @title Multicore BART
#' 
#' @description This function performs Bayesian Additive Regression Trees (BART) using multiple cores.
#' 
#' @param x.train Training data matrix.
#' @param y.train Response variable for training data.
#' @param x.test Test data matrix (optional).
#' @param sparse Logical. If TRUE, sparse matrix is used.
#' @param theta Parameter for the model.
#' @param omega Parameter for the model.
#' @param a Parameter for the model.
#' @param b Parameter for the model.
#' @param augment Logical. If TRUE, augmentation is used.
#' @param rho Parameter for the model.
#' @param xinfo Information about the model matrix.
#' @param usequants Logical. If TRUE, quantiles are used for cut points.
#' @param cont Logical. If TRUE, continuous variables are treated as such.
#' @param rm.const Logical. If TRUE, constant variables are removed.
#' @param sigest Initial estimate of the error variance.
#' @param sigdf Degrees of freedom for the error variance.
#' @param sigquant Quantile for the error variance.
#' @param k Parameter for the model.
#' @param power Parameter for the model.
#' @param base Parameter for the model.
#' @param sigmaf Parameter for the model.
#' @param lambda Parameter for the model.
#' @param fmean Mean of the response variable.
#' @param w Weights for the response variable.
#' @param ntree Number of trees to grow.
#' @param numcut Number of cut points for each variable.
#' @param ndpost Number of posterior samples to draw.
#' @param nskip Number of samples to skip.
#' @param keepevery Keep every nth sample.
#' @param printevery Print progress every nth sample.
#' @param keeptrainfits Logical. If TRUE, training fits are kept.
#' @param transposed Logical. If TRUE, the data is transposed.
#' @param mc.cores Number of cores to use for parallel processing.
#' @param nice Nice level for the process.
#' @param seed Random seed for reproducibility.
#' 
#' @return A list containing the posterior samples, predictions, and other information.
#' 
#' @export
#' 
#' @importFrom parallel detectCores
#' @importFrom tools psnice
#' @importFrom abind abind
#' 
mc.mebart <- function(
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
                mebart(
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



    if (mc.cores == 1 | attr(post, "class") != "mebart") {
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

        attr(post, "class") <- "mebart"

        return(post)
    }
}
