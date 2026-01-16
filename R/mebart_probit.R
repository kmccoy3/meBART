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




#' @title Fit a meBART Model with \emph{Binary} Output Data
#'
#' @description This function fits a Bayesian Additive Regression Tree (BART) model which also also models measurement error in the explanatory variables. 
#' This function specifically handles \emph{binary} output data.
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
#' \item{x_draws}{Tensor with posterior draws for the training data.}
#' \item{acceptances}{Whether or not the proposal was accepted.}
#' \item{proc.time}{The time taken for the C++ code to run.}
#' \item{varcount.mean}{The mean of the variable counts.}
#' \item{varprob.mean}{The mean of the variable probabilities.}
#' \item{rm.const}{A vector specifying which columns were removed from `x.train`.}
#' @export
#' @importFrom stats lm
#'
#' @examples
#' # Example usage:
#' x.train <- matrix(rnorm(1000), ncol = 10)
#' y.train <- rbinom(100, 1, 0.5)
#' x.test <- matrix(rnorm(100), ncol = 10)
#' y.test <- rbinom(10, 1, 0.5)
#' mdl <- mebart_probit(x.train, y.train, x.test,
#'                ndpost = 500,
#'                ntree = 50,
#'                meas_error_sigma = diag(10),
#'                x_mu = matrix(0, nrow=10, ncol=1),
#'                x_sigma = diag(10))
#'  
#' preds <- mdl$prob.test.mean
#' print(paste0("Test Accuracy: ", mean((preds > 0.5) == y.test)))
#'
mebart_probit = function(x.train,
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
                 k = 2.0,
                 power = 2.0,
                 base = .95,
                 binaryOffset = NULL,
                 ntree = 50L,
                 numcut = 100L,
                 ndpost = 1000L,
                 nskip = 100L,
                 keepevery = 1L,
                 nkeeptrain = ndpost,
                 nkeeptest = ndpost,
                 ##nkeeptestmean=ndpost,
                 nkeeptreedraws = ndpost,
                 printevery = 100L,
                 transposed = FALSE, # used if called by mc.wbart
                 proposal_sigma = meas_error_sigma, # standard deviation of the proposal distribution
                 meas_error_sigma, # standard deviation of the measurement error
                 x_mu,
                 x_sigma
                 ) {
    #--------------------------------------------------
    #data
    n = length(y.train)
    
    if (length(binaryOffset) == 0)
        binaryOffset = qnorm(mean(y.train))
    
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
        ##         x.test = t(bartModelMatrix(x.test[ , temp$rm.const]))
        if (length(x.test) > 0) {
            x.test = bartModelMatrix(x.test)
            x.test = t(x.test[, temp$rm.const])
        }
        rm.const <- temp$rm.const
        grp <- temp$grp
        rm(temp)
    }
    else {
        rm.const <- NULL
        grp <- NULL
    }
    
    if (n != ncol(x.train))
        stop('The length of y.train and the number of rows in x.train must be identical')
    
    p = nrow(x.train)
    np = ncol(x.test)
    if (length(rho) == 0)
        rho <- p
    if (length(rm.const) == 0)
        rm.const <- 1:p
    if (length(grp) == 0)
        grp <- 1:p


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
    
    #--------------------------------------------------
    #set  nkeeps for thinning
    if ((nkeeptrain != 0) &
        ((ndpost %% nkeeptrain) != 0)) {
        nkeeptrain = ndpost
        cat('*****nkeeptrain set to ndpost\n')
    }
    if ((nkeeptest != 0) &
        ((ndpost %% nkeeptest) != 0)) {
        nkeeptest = ndpost
        cat('*****nkeeptest set to ndpost\n')
    }
    ## if((nkeeptestmean!=0) & ((ndpost %% nkeeptestmean) != 0)) {
    ##    nkeeptestmean=ndpost
    ##    cat('*****nkeeptestmean set to ndpost\n')
    ## }
    if ((nkeeptreedraws != 0) &
        ((ndpost %% nkeeptreedraws) != 0)) {
        nkeeptreedraws = ndpost
        cat('*****nkeeptreedraws set to ndpost\n')
    }
    #--------------------------------------------------
    #prior
    ## nu=sigdf
    ## if(is.na(lambda)) {
    ##    if(is.na(sigest)) {
    ##       if(p < n) {
    ##          df = data.frame(t(x.train),y.train)
    ##          lmf = lm(y.train~.,df)
    ##          rm(df)
    ##          sigest = summary(lmf)$sigma
    ##       } else {
    ##          sigest = sd(y.train)
    ##       }
    ##    }
    ##    qchi = qchisq(1.0-sigquant,nu)
    ##    lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior
    ## }
    
    ## if(is.na(sigmaf)) {
    ##    tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree));
    ## } else {
    ##    tau = sigmaf/sqrt(ntree)
    ## }
    #--------------------------------------------------
    ptm <- proc.time()
    #call
    res = .Call(
        "cpbart",
        n,
        #number of observations in training data
        p,
        #dimension of x
        np,
        #number of observations in test data
        x.train,
        #p*n training data x
        y.train,
        #n*1 training data y
        x.test,
        #p*np test data x
        ntree,
        numcut,
        ndpost * keepevery,
        nskip,
        power,
        base,
        binaryOffset,
        3 / (k * sqrt(ntree)),
        sparse,
        theta,
        omega,
        grp,
        a,
        b,
        rho,
        augment,
        nkeeptrain,
        nkeeptest,
        ##nkeeptestmean,
        nkeeptreedraws,
        printevery,
        ##treesaslists,
        xinfo, # cutpoints, now specified
        proposal_sigma, # standard deviation of the proposal distribution
        meas_error_sigma, # standard deviation of the measurement error
        x_mu,
        x_sigma
    )
    
    res$proc.time <- proc.time() - ptm
    
    if (nkeeptrain > 0) {
        ##res$yhat.train.mean <- NULL
        ##res$yhat.train.mean = res$yhat.train.mean+binaryOffset
        res$yhat.train = res$yhat.train + binaryOffset
        res$prob.train = pnorm(res$yhat.train)
        res$prob.train.mean <- apply(res$prob.train, 2, mean)
    } else {
        res$yhat.train <- NULL
        ##res$yhat.train.mean <- NULL
    }
    
    if (np > 0) {
        ##res$yhat.test.mean <- NULL
        ##res$yhat.test.mean = res$yhat.test.mean+binaryOffset
        res$yhat.test = res$yhat.test + binaryOffset
        res$prob.test = pnorm(res$yhat.test)
        res$prob.test.mean <- apply(res$prob.test, 2, mean)
    } else {
        res$yhat.test <- NULL
        ##res$yhat.test.mean <- NULL
    }
    
    if (nkeeptreedraws > 0)
        ## & !treesaslists)
        names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
    
    dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
    dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
    res$varcount.mean <- apply(res$varcount, 2, mean)
    res$varprob.mean <- apply(res$varprob, 2, mean)
    res$rm.const <- rm.const
    res$binaryOffset = binaryOffset
    attr(res, 'class') <- 'mebart_probit'
    return(res)
}