

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

wbart = function(x.train, # explanatory variables for training data, matrix or dataframe
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
                 transposed = FALSE) # used if called by mc.wbart
{
    #--------------------------------------------------
    #data
    n = length(y.train)
    
    if (!transposed) { # if being called wbart, skipped if being called by mc.wbart
        temp = bartModelMatrix(
            x.train,
            numcut,
            usequants = usequants, # if false makes cutpoints uniform, if true uniform quantiles are used
            cont = cont, # whether or not to assume all variables are continuous
            xinfo = xinfo, # cutpoints, if you want to specify them
            rm.const = rm.const # remove constant columns from x.train
        )
        x.train = t(temp$X) # transpose data matrix, has been converted to matrix
        numcut = temp$numcut # number of cutpoints, vector of each variable
        xinfo = temp$xinfo # cutpoints, now specified
        if (length(x.test) > 0) {
            x.test = bartModelMatrix(x.test)
            x.test = t(x.test[, temp$rm.const])
        }
        rm.const <- temp$rm.const # removed variables or all if kept
        grp <- temp$grp # counts out which variables are unique / which ones are factors
        rm(temp) # clears temp from R's workspace 
    }
    else {
        rm.const <- NULL    
        grp <- NULL
    }
    
    if (n != ncol(x.train)) # check if the number of observations in x.train is equal to the length of y.train, remember x.train is transposed
        stop('The length of y.train and the number of rows in x.train must be identical')
    
    p = nrow(x.train)
    np = ncol(x.test)
    if (length(rho) == 0)
        rho = p
    if (length(rm.const) == 0)
        rm.const <- 1:p
    if (length(grp) == 0)
        grp <- 1:p
    
    ##if(p>1 & length(numcut)==1) numcut=rep(numcut, p)
    
    y.train = y.train - fmean # centers output data
    #------------------------------------------------------------------------------------------------------------
    #set nkeeps for thinning ######### Basically makes sure these variables dont conflict with one another
    if ((nkeeptrain != 0) & ((ndpost %% nkeeptrain) != 0)) { # if 
        nkeeptrain = ndpost
        cat('*****nkeeptrain set to ndpost\n') # cat is just a simpler version of print
    }
    if ((nkeeptest != 0) & ((ndpost %% nkeeptest) != 0)) {
        nkeeptest = ndpost
        cat('*****nkeeptest set to ndpost\n')
    }
    if ((nkeeptestmean != 0) & ((ndpost %% nkeeptestmean) != 0)) {
        nkeeptestmean = ndpost
        cat('*****nkeeptestmean set to ndpost\n')
    }
    if ((nkeeptreedraws != 0) & ((ndpost %% nkeeptreedraws) != 0)) {
        nkeeptreedraws = ndpost
        cat('*****nkeeptreedraws set to ndpost\n')
    }
    #----------------------------------------------------------------------------------------------
    #prior
    nu = sigdf # nu / DOF in sigma prior
    if (is.na(lambda)) { # NA by default
        if (is.na(sigest)) { # NA by default
            if (p < n) {
                df = data.frame(t(x.train), y.train)
                lmf = lm(y.train ~ ., df) # trains linear model
                sigest = summary(lmf)$sigma # gets sigma from linear model
            } else {
                sigest = sd(y.train)
            }
        }
        qchi = qchisq(1.0 - sigquant, nu) # gets quantile from chi squared distribution
        lambda = (sigest * sigest * qchi) / nu # lambda parameter for sigma prior
    } else {
        sigest = sqrt(lambda)
    }
    
    if (is.na(sigmaf)) { # NA by default
        tau = (max(y.train) - min(y.train)) / (2 * k * sqrt(ntree)) 
    } else {
        tau = sigmaf / sqrt(ntree)
    }
    #-------------------------------------------------------------------------------------------------------
    ptm <- proc.time() # how much real and CPU time (in seconds) the currently running R process has already taken.
    #call
    res = .Call(
        "cwbart",
        n, #number of observations in training data
        p, #dimension of x
        np, #number of observations in test data
        x.train, #pxn training data x
        y.train, #pxn training data x
        x.test,  #p*np test data x
        ntree, #number of trees
        numcut, #number of cutpoints
        ndpost * keepevery, # total number of draws
        nskip, # number of observations for burn-in
        power, # tree depth prior, beta in CGM
        base, # tree depth prior, alpha in CGM
        tau,  # sigma from mu prior, now scaled by number of trees
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
        xinfo # cutpoints, now specified
    )
    
    res$proc.time <- proc.time() - ptm # how much time C++ code has taken
    
    # Centering variable, add back in
    res$mu = fmean
    res$yhat.train.mean = res$yhat.train.mean + fmean
    res$yhat.train = res$yhat.train + fmean
    res$yhat.test.mean = res$yhat.test.mean + fmean
    res$yhat.test = res$yhat.test + fmean



    if (nkeeptreedraws > 0)
        names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
    dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
    dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
    ##res$nkeeptreedraws=nkeeptreedraws
    res$varcount.mean <- apply(res$varcount, 2, mean)
    res$varprob.mean <- apply(res$varprob, 2, mean)
    res$rm.const <- rm.const
    attr(res, 'class') <- 'wbart'
    return(res)
}
