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


#' @title Fit a meBART Model with \emph{Binary} Output Data (Multi-core Version)
#'
#' @description This function fits a Bayesian Additive Regression Tree (BART) model which also also models measurement error in the explanatory variables. This function specifically handles \emph{binary} output data.
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
#' @param k Prior parameter for the mean of the output data.
#' @param power Power parameter for the BART tree-depth prior (beta in original BART paper).
#' @param base Base parameter for the BART tree-depth prior (alpha in original BART paper).
#' @param binaryOffset Offset for binary outcomes. Defaults to Phi(mean(y.train)).
#' @param ntree The number of trees to use in the model.
#' @param numcut The number of cutpoints to consider for each variable.
#' @param ndpost The number of posterior draws to return.
#' @param nskip The number of observations for burn-in (MCMC).
#' @param keepevery The thinning parameter for MCMC.
#' @param nkeeptrain Number of training predictions to keep.
#' @param nkeeptest Number of test predictions to keep.
#' @param nkeeptreedraws Number of tree draws to keep.
#' @param printevery How often to print progress during MCMC.
#' @param transposed Is the input data transposed? Used if called by `mc.wbart`.
#' @param mc.cores The number of cores to use for parallel processing.
#' @param nice The niceness level for the parallel processes (default is 19, which is the lowest priority).
#' @param seed The random seed for reproducibility.
#' @param proposal_sigma A numeric value (default is `meas_error_sigma`). Covariance matrix of the proposal distribution.
#' @param meas_error_sigma A numeric matrix. Covariance matrix of the measurement error, must be square and match the number of columns in `x.train`.
#' @param x_mu A numeric column vector. Mean of the explanatory variables.
#' @param x_sigma A numeric matrix. Covariance matrix of the explanatory variables.
#'
#' @return A list with the following elements:
#' \item{yhat.train}{An (ndpost * n) matrix of the posterior predictions for the training data.}
#' \item{yhat.test}{An (ndpost * np) matrix of the posterior predictions for the test data.}
#' \item{varcount}{An (ndpost * p) matrix of variable counts in the model.}
#' \item{varprob}{An (ndpost * p) matrix of variable probabilities in the model.}
#' \item{treedraws}{A list of posterior tree draws.}
#' \item{x_draws}{Tensor with posterior draws for the training data.}
#' \item{acceptances}{Whether or not the proposal was accepted.}
#' \item{proc.time}{The time taken for the C++ code to run.}
#' \item{prob.train}{An (ndpost * n) matrix of the posterior probabilities for the training data.}
#' \item{prob.train.mean}{The mean of the posterior probabilities for the training data.}
#' \item{prob.test}{An (ndpost * np) matrix of the posterior probabilities for the test data.}
#' \item{prob.test.mean}{The mean of the posterior probabilities for the test data.}
#' \item{varcount.mean}{The mean of the variable counts.}
#' \item{varprob.mean}{The mean of the variable probabilities.}
#' \item{rm.const}{A vector specifying which columns were removed from `x.train`.}
#' \item{binaryOffset}{The offset used for binary outcomes.}
#' @export
#' @importFrom stats lm
#'
#' @examples
#' # Example usage:
#' x.train <- matrix(rnorm(1000), ncol = 10)
#' y.train <- rbinom(100, 1, 0.5)
#' x.test <- matrix(rnorm(100), ncol = 10)
#' y.test <- rbinom(10, 1, 0.5)
#' mdl <- mc.mebart_probit(x.train, y.train, x.test,
#'                ndpost = 500,
#'                ntree = 50,
#'                meas_error_sigma = diag(10),
#'                x_mu = matrix(0, nrow=10, ncol=1),
#'                x_sigma = diag(10),
#'                mc.cores = 2)
#'  
#' preds <- mdl$prob.test.mean
#' print(paste0("Test Accuracy: ", mean((preds > 0.5) == y.test)))
#' 
mc.mebart_probit <- function(x.train,
                     y.train,
                     x.test = matrix(0.0, 0L, 0L),
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
                     k = 2.0,
                     ## BEWARE: do NOT use k for other purposes below
                     power = 2.0,
                     base = 0.95,
                     binaryOffset = NULL,
                     ntree = 50L,
                     numcut = 100L,
                     ndpost = 1000L,
                     nskip = 100L,
                     keepevery = 1L,
                     printevery = 100L,
                     keeptrainfits = TRUE,
                     transposed = FALSE,
                     ##    treesaslists=FALSE,
                     mc.cores = 2L,
                     nice = 19L,
                     seed = 99L,
                     proposal_sigma = meas_error_sigma, # standard deviation of the proposal distribution
                     meas_error_sigma, # standard deviation of the measurement error
                     x_mu,
                     x_sigma
                     )
{
    if (.Platform$OS.type != 'unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')
    
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()
    
    if (!transposed) {
        temp = bartModelMatrix(
            x.train,
            numcut,
            usequants = usequants,
            cont = cont,
            xinfo = xinfo,
            rm.const = rm.const
        )
        x.train = t(temp$X)
        numcut = temp$numcut
        xinfo = temp$xinfo
        ## if(length(x.test)>0)
        ##     x.test = t(bartModelMatrix(x.test[ , temp$rm.const]))
        if (length(x.test) > 0) {
            x.test = bartModelMatrix(x.test)
            x.test = t(x.test[, temp$rm.const])
        }
        rm.const <- temp$rm.const
        rm(temp)
    }
    
    mc.cores.detected <- detectCores()
    
    if (mc.cores > mc.cores.detected)
        mc.cores <- mc.cores.detected
    ## warning(paste0('The number of cores requested, mc.cores=', mc.cores,
    ##                ',\n exceeds the number of cores detected via detectCores() ',
    ##                'which yields ', mc.cores.detected, ' .'))
    
    mc.ndpost <- ceiling(ndpost / mc.cores)
    ## mc.ndpost <- ((ndpost %/% mc.cores) %/% keepevery)*keepevery
    
    ## while(mc.ndpost*mc.cores<ndpost) mc.ndpost <- mc.ndpost+keepevery
    
    ## mc.nkeep <- mc.ndpost %/% keepevery
    
    for (i in 1:mc.cores) {
        parallel::mcparallel({
            psnice(value = nice)
            
            mebart_probit(
                x.train = x.train,
                y.train = y.train,
                x.test = x.test,
                sparse = sparse,
                theta = theta,
                omega = omega,
                a = a,
                b = b,
                augment = augment,
                rho = rho,
                xinfo = xinfo,
                k = k,
                power = power,
                base = base,
                binaryOffset = binaryOffset,
                ntree = ntree,
                numcut = numcut,
                ndpost = mc.ndpost,
                nskip = nskip,
                keepevery = keepevery,
                ## nkeeptrain=mc.nkeep, nkeeptest=mc.nkeep,
                ## nkeeptestmean=mc.nkeep, nkeeptreedraws=mc.nkeep,
                printevery = printevery,
                transposed = TRUE,
                    proposal_sigma = meas_error_sigma, # standard deviation of the proposal distribution
                   meas_error_sigma = meas_error_sigma, # standard deviation of the measurement error
                   x_mu = x_mu,
                   x_sigma = x_sigma
            ) 
        }, ##treesaslists=treesaslists)},
        silent = (i != 1))
        ## to avoid duplication of output
        ## capture stdout from first posterior only
    }
    
    post.list <- parallel::mccollect()
    
    post <- post.list[[1]]
    
    if (mc.cores == 1 | attr(post, 'class') != 'mebart_probit')
        return(post)
    else {
        if (class(rm.const)[1] != 'logical')
            post$rm.const <- rm.const
        
        post$ndpost <- mc.cores * mc.ndpost
        
        p <- nrow(x.train[post$rm.const, ])
        ##p <- nrow(x.train[ , post$rm.const])
        
        ## if(length(rm.const)==0) rm.const <- 1:p
        ## post$rm.const <- rm.const
        
        old.text <- paste0(as.character(mc.ndpost),
                           ' ',
                           as.character(ntree),
                           ' ',
                           as.character(p))
        old.stop <- nchar(old.text)
        
        post$treedraws$trees <- sub(
            old.text,
            paste0(
                as.character(post$ndpost),
                ' ',
                as.character(ntree),
                ' ',
                as.character(p)
            ),
            post$treedraws$trees
        )
        
        keeptestfits <- length(x.test) > 0
        
        for (i in 2:mc.cores) {
            if (keeptrainfits) {
                post$yhat.train <- rbind(post$yhat.train, post.list[[i]]$yhat.train)
                post$prob.train <- rbind(post$prob.train, post.list[[i]]$prob.train)
            }
            
            if (keeptestfits) {
                post$yhat.test <- rbind(post$yhat.test, post.list[[i]]$yhat.test)
                post$prob.test <- rbind(post$prob.test, post.list[[i]]$prob.test)
            }
            
            post$varcount <- rbind(post$varcount, post.list[[i]]$varcount)
            post$varprob <- rbind(post$varprob, post.list[[i]]$varprob)
            
            post$treedraws$trees <- paste0(
                post$treedraws$trees,
                substr(
                    post.list[[i]]$treedraws$trees,
                    old.stop + 2,
                    nchar(post.list[[i]]$treedraws$trees)
                )
            )
            
            ## if(treesaslists) post$treedraws$lists <-
            ##                      c(post$treedraws$lists, post.list[[i]]$treedraws$lists)
        }
        
        ## if(length(post$yhat.train.mean)>0)
        ##     post$yhat.train.mean <- apply(post$yhat.train, 2, mean)
        
        ## if(length(post$yhat.test.mean)>0)
        ##     post$yhat.test.mean <- apply(post$yhat.test, 2, mean)
        
        if (length(post$prob.train.mean) > 0)
            post$prob.train.mean <- apply(post$prob.train, 2, mean)
        
        if (length(post$prob.test.mean) > 0)
            post$prob.test.mean <- apply(post$prob.test, 2, mean)
        
        post$varcount.mean <- apply(post$varcount, 2, mean)
        post$varprob.mean <- apply(post$varprob, 2, mean)
        
        attr(post, 'class') <- 'mebart_probit'
        
        return(post)
    }
}