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


#' @title Single-core + Multi-core (OpenMP) meBART Predictions
#'
#' @description This function performs single-core or multi-core (OpenMP) predictions for meBART 
#' models. This function is also called by `mc.pmebart.R`.
#'
#' @param x.test A matrix of test data.
#' @param treedraws A list containing the tree draws. (Part of the meBART model)
#' @param mu Offset for the predictions.
#' @param mc.cores Number of cores to use for parallel processing.
#' @param transposed Logical. Is x.test transposed?
#' @param dodraws Logical. If TRUE, all posterior draws are returned; if FALSE, only the mean is returned.
#' @param nice Nice level for the process. I'm pretty sure this is not used in the code.
#'
#' @return Prediction matrix of size (ndpost x n).
#'
#' @keywords internal

pmebart <- function(x.test,
                    treedraws,
                    mu = 0,
                    mc.cores = 1L,
                    transposed = FALSE,
                    dodraws = TRUE,
                    nice = 19L
                    ) {

    # Transpose x.test if it hasn't been already 
    if (!transposed){
        x.test <- t(bartModelMatrix(x.test))
    }
    
    # Check that new data conforms to the trained model
    p <- length(treedraws$cutpoints)
    if (p != nrow(x.test)) {
        stop(paste0("The number of columns in x.test must be equal to ", p))
    }
    
    # Call the C++ function for predictions
    res <- .Call(
        "cpmebart",
        treedraws,
        x.test,
        mc.cores
    )

    # Return all posterior draws or just the means
    if (dodraws) {
        return(res$yhat.test + mu)
    } else {
        return(apply(res$yhat.test, 2, mean) + mu)
    }

}
