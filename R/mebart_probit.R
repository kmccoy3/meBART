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
#' @importFrom stats lm pnorm qnorm
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