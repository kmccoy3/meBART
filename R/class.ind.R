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


#' @title Class Indicator Matrix
#' 
#' @description This function one hot encodes a factor to create a class indicator matrix.
#' 
#' @param cl A factor or vector of class labels.
#' 
#' @return A (n x num_levels) matrix with rows corresponding to observations and columns corresponding to each unique level of the factor.
#' 
#' @keywords internal

class.ind <- function(cl)
{
    # Get number of observations
    n <- length(cl)

    # If not already a factor, convert to factor
    cl <- as.factor(cl)

    # Create zero-matrix with n rows and number of levels in cl columns
    x <- matrix(0, n, length(levels(cl)) )

    # Set the appropriate entries to 1
    x[(1L:n) + n*(unclass(cl)-1L)] <- 1

    # Set row and column names, rownames are usually NULL
    dimnames(x) <- list(names(cl), levels(cl))

    # Return the class indicator matrix
    x
}