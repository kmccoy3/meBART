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


#' @title Multicore PROBIT BART
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
mc.pbart <- function(x.train,
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
                     seed = 99L)
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
            
            pbart(
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
                transposed = TRUE
            )
        }, ##treesaslists=treesaslists)},
        silent = (i != 1))
        ## to avoid duplication of output
        ## capture stdout from first posterior only
    }
    
    post.list <- parallel::mccollect()
    
    post <- post.list[[1]]
    
    if (mc.cores == 1 | attr(post, 'class') != 'pbart')
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
        
        attr(post, 'class') <- 'pbart'
        
        return(post)
    }
}