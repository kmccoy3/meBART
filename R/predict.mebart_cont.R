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


#' @title Make Predictions for meBART Model with Continuous Outcomes
#'
#' @description This function makes predictions using a meBART model with \emph{continuous} outcomes.
#'
#' @param object A meBART object.
#' @param newdata A matrix of new data for prediction.
#' @param mc.cores Number of cores to use for parallel processing.
#' @param openmp Logical. If TRUE, OpenMP is used for parallel processing.
#' @param ... Additional arguments.
#'
#' @return A matrix of predictions with size (ndpost x n)
#'
#' @export
#' @importFrom parallel detectCores
#'
#' @examples
#' # Example usage:
#' set.seed(0)
#' x.train <- matrix(rnorm(1000), ncol = 10)
#' y.train <- rnorm(100)
#' mdl <- meBART::mebart(x.train, y.train,
#'               ndpost = 500,
#'               ntree = 50,
#'               meas_error_sigma = diag(10),
#'               x_mu = matrix(0, nrow=10, ncol=1),
#'               x_sigma = diag(10))
#' preds <- predict(mdl, newdata = x.train)

predict.mebart_cont <- function(object, # A meBART object.
                           newdata, # A matrix of new data for prediction.
                           mc.cores = 1, # Number of cores to use for parallel processing.
                           openmp = (mc.cores.openmp() > 0), # If TRUE, OpenMP is used for parallel processing.
                           ...) {
    
    # Check that new data conforms to the trained model
    p <- length(object$treedraws$cutpoints)
    if (p != ncol(newdata)) {
        stop(paste0("The number of columns in newdata must be equal to ", p))
    }
    
    # Detect the number of cores available
    if (.Platform$OS.type == "unix") {
        mc.cores.detected <- parallel::detectCores()
    } else {
        mc.cores.detected <- NA
    }
    
    # Will max out cores at the detected number
    if (!is.na(mc.cores.detected) && mc.cores > mc.cores.detected) {
        mc.cores <- mc.cores.detected
    }
    
    # Either use single core predictions or multicore predictions in OpenMP / C++
    if (.Platform$OS.type != "unix" || openmp || mc.cores == 1) {
        call <- pmebart
    } else { 
        call <- mc.pmebart # Multicore predictions using parallel package in R
    }
    
    # Not sure if this is necessary, but it is in the original code
    if (length(object$mu) == 0) {
        object$mu <- object$offset
    }
    
    # Call the prediction function with the new data and model parameters
    return(call(
        newdata,
        object$treedraws,
        mc.cores = mc.cores,
        mu = object$mu,
        ...
    ))
}
