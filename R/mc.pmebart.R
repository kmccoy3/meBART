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


#' @title Multi-core (parallel) meBART Predictions
#'
#' @description This function performs multi-core (parallel) predictions for meBART models.
#'
#' @param x.test A matrix of test data.
#' @param treedraws A list containing the tree draws. (Part of the meBART model)
#' @param mu Offset for the predictions.
#' @param mc.cores Number of cores to use for parallel processing.
#' @param transposed Logical. Is x.test transposed?
#' @param dodraws Logical. If TRUE, all posterior draws are returned; if FALSE, only the mean is returned.
#' @param nice Nice level for the process.
#'
#' @return Prediction matrix of size (ndpost x n).
#'
#' @keywords internal
#' @importFrom parallel detectCores
#' @importFrom tools psnice

mc.pmebart <- function(x.test,
                       treedraws,
                       mu = 0,
                       mc.cores = 2L,
                       transposed = FALSE,
                       dodraws = TRUE,
                       nice = 19L
                       ) {

    # Check if parallel processing is available
    if (.Platform$OS.type != "unix") {
        stop("parallel::mcparallel/mccollect do not exist on windows")
    }
    
    # Transpose x.test if it hasn't been already
    if (!transposed) {
        x.test <- t(bartModelMatrix(x.test))
    }
    
    # Check that new data conforms to the trained model
    p <- length(treedraws$cutpoints)
    if (p != nrow(x.test)) {
        stop(paste0("The number of columns in x.test must be equal to ", p))
    }
    
    # Detect the number of cores available
    mc.cores.detected <- parallel::detectCores()
    
    # Will max out cores at the detected number
    if (!is.na(mc.cores.detected) && mc.cores > mc.cores.detected) {
        mc.cores <- mc.cores.detected
    }
    
    # Get number of test observations, limit mc.cores to that number
    K <- ncol(x.test)
    if (K < mc.cores) mc.cores <- K

    # How many observations to predict per core
    k <- K %/% mc.cores - 1

    # Initialize j, end index
    j <- K
    
    # Loop through cores
    for (i in 1:mc.cores) {

        # Update starting index for this core
        if (i == mc.cores) {
            h <- 1 # Start index for the last core
        } else {
            h <- j - k # Ending index minus k
        }

        # Issue a parallel task for this core
        parallel::mcparallel({
            tools::psnice(value = nice);
            pmebart(matrix(x.test[, h:j], nrow = p, ncol = j - h + 1),
                    treedraws,
                    mu,
                    1,
                    TRUE)
        }, silent = (i != 1))

        j <- h - 1 # new end index for next core
    }
    
    # Get predictions from all cores
    pred.list <- parallel::mccollect()
    pred <- pred.list[[1]] # First core's prediction

    # Check the type of prediction returned
    type <- class(pred)[1]
    if (type == "list") {
        pred <- pred[[1]]
    } else if (type != "matrix") {
        return(pred.list)
    } ## likely error messages
    
    # Combine predictions from all cores
    if (mc.cores > 1) {
        for (i in 2:mc.cores) {
            if (type == "list") {
                pred <- cbind(pred, pred.list[[i]][[1]])
            } else {
                pred <- cbind(pred, pred.list[[i]])
            }
        }
    }
    
    # Return all posterior draws or just the means
    if (dodraws) {
        return(pred) # Recall mu has already been added in pmebart
    } else {
        return(apply(pred, 2, mean))
    }

}
