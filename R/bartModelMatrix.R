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


#' @title Create a matrix out of a vector or data.frame
#' 
#' @description The external BART functions operate on matrices in memory. Therefore, if the user 
#' submits a vector or data.frame, then this function converts it to a matrix. Also, it determines 
#' the number of cutpoints necessary and their values for each column when asked to do so.
#' 
#' @param X A vector or data.frame to create the matrix from.
#' @param numcut Integer. The maximum number of cutpoints to consider. If numcut=0, then just 
#' return a matrix; otherwise, return a list containing a matrix X, a vector numcut and a list xinfo.
#' @param usequants Logical. If usequants is FALSE, then the cutpoints in xinfo are generated 
#' uniformly; otherwise, if TRUE, then quantiles are used for the cutpoints.
#' @param type Determines which quantile algorithm is employed.
#' @param rm.const Logical. Whether or not to remove constant variables.
#' @param cont Logical. Whether or not to assume all variables are continuous.
#' @param xinfo You can provide the cutpoints to BART or let BART choose them for you. To provide 
#' them, use the xinfo argument to specify a list (matrix) where the items (rows) are the 
#' covariates and the contents of the items (columns) are the cutpoints.
#' 
#' @return A list containing the cleaned input data, number of cut points, removed variables, and cut point information.
#' 
#' @keywords internal

bartModelMatrix <- function(X,
                            numcut = 0L,
                            usequants = FALSE,
                            type = 7,
                            rm.const = FALSE,
                            cont = FALSE,
                            xinfo = NULL
                            ) {

    # Data.frame or matrix, etc.
    X.class <- class(X)[1]

    # If single column of factor
    if (X.class == "factor") {
        X.class <- "data.frame"
        X <- data.frame(X = X)
    }

    # Initialize variable groupings
    grp <- NULL

    # ==============================================================================================
    # Clean input data, coerce to matrix
    # ==============================================================================================

    if (X.class == "data.frame") {

        # Number of columns
        p <- dim(X)[2]
        # Column names
        xnm <- names(X)

        # Loop over features
        for (i in 1:p) {
            if (is.factor(X[[i]])) {
                # Xtemp <- class.ind(X[[i]])
                # colnames(Xtemp) <- paste(xnm[i], 1:ncol(Xtemp), sep = "")
                # X[[i]] <- Xtemp
                # grp <- c(grp, rep(i, ncol(Xtemp)))
                stop("Factor inputs are not currently supported for meBART. All features must be numeric.")
            } else {
                X[[i]] <- cbind(X[[i]]) # Clean each column
                colnames(X[[i]]) <- xnm[i]
                grp <- c(grp, i) # Variable group (used for one-hot encoding)
            }
        }
        # Clean X
        Xtemp <- cbind(X[[1]])
        if (p > 1) for (i in 2:p) Xtemp <- cbind(Xtemp, X[[i]])
        X <- Xtemp
    } else if (X.class == "numeric" | X.class == "integer") {
        # If a single numeric vector, convert to matrix
        X <- cbind(as.numeric(X))
        grp <- 1
    } else if (X.class == "NULL") {
        return(X)
    } else if (X.class != "matrix") {
        stop("Expecting either a factor, a vector, a matrix or a data.frame")
    }

    # ==============================================================================================
    # Calculate cutpoints
    # ==============================================================================================

    # Get dimensions of data
    N <- nrow(X)
    p <- ncol(X)

    # Initialize variables
    xinfo. <- matrix(nrow = p, ncol = numcut)
    nc <- numcut
    rm.vars <- c()

    if (N > 0 & p > 0 & (rm.const | numcut[1] > 0)) {
        for (j in 1:p) {

            # Get data type for each column
            X.class <- class(X[1, j])[1]

            if (X.class == "numeric" | X.class == "integer") {

                # Get sorted, unique values (eligible cut points)
                xs <- unique(sort(X[, j]))
                k <- length(xs)
                nc[j] <- numcut

                if (k %in% 0:1) {
                    # If only one unique value, remove that column
                    rm.vars <- c(rm.vars, -j)
                    nc[j] <- 1
                    if (k == 0) xs <- NA
                } else if (cont) {
                    xs <- seq(xs[1], xs[k], length.out = numcut + 2)[-c(1, numcut + 2)]
                } else if (k < numcut) {
                    xs <- 0.5 * (xs[1:(k - 1)] + xs[2:k])
                    nc[j] <- k - 1
                } else if (usequants) {
                    xs <- quantile(X[, j], type = type, probs = (0:(numcut + 1)) / (numcut + 1))[-c(1, numcut + 2)]
                    names(xs) <- NULL
                } else {
                    xs <- seq(xs[1], xs[k], length.out = numcut + 2)[-c(1, numcut + 2)]
                }
            } else {
                stop(paste0("Variables of type ", X.class, " are not supported"))
            }

            # Fill in xinfo with cut points
            xinfo.[j, 1:nc[j]] <- xs
        }
    }

    # ==============================================================================================
    # Clean Xinfo and remove constant variables
    # ==============================================================================================

    # If xinfo is provided, clean it and convert it to a matrix
    if (length(xinfo) > 0) {
        if (is.list(xinfo)) {
            for (j in 1:p) {
                # If xinfo is a list, convert it to a matrix
                xinfo.[j, 1:length(xinfo[[j]])] <- xinfo[[j]]
            }
        } else if (is.matrix(xinfo)) {
            xinfo. <- xinfo
        } else {
            stop("Only a list or a matrix can be provided for xinfo")
        }

        for (j in 1:p) {
            # Number of cut points for each variable
            nc[j] <- sum(!is.na(xinfo.[j, ]))
        }
    }

    xinfo <- xinfo.

    # Convert data.frame to matrix
    X <- data.matrix(X)

    if (rm.const & length(rm.vars) > 0) {
        # Remove constant variables
        X <- X[, rm.vars]
        nc <- nc[rm.vars]
        xinfo <- xinfo[rm.vars, ]
    } else if (length(rm.vars) == 0) {
        # If no variables were removed, set rm.vars to 1:p (i.e., all variables used)
        rm.vars <- 1:p
    }

    # Rows of xinfo are the variables, columns are the cut points
    dimnames(xinfo) <- list(dimnames(X)[[2]], NULL)

    # ==============================================================================================
    # Return the model matrix and other information
    # ==============================================================================================

    # Return just X or full list, depending on where this function is called
    if (numcut == 0) {
        return(X)
    } else {
        return(list(
            X = X, 
            numcut = as.integer(nc), 
            rm.const = rm.vars,
            xinfo = xinfo, 
            grp = grp
        ))
    }
}
