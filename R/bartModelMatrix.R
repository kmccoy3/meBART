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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/GPL-2


#' @title BART Model Matrix
#' @description This function creates a model matrix for BART.
#' 
#' @param X A data frame, matrix, or vector.
#' @param numcut The number of cut points for each variable.
#' @param usequants Logical. If TRUE, quantiles are used for cut points.
#' @param type The type of quantile to use.
#' @param rm.const Logical. If TRUE, constant variables are removed.
#' @param cont Logical. If TRUE, continuous variables are treated as such.
#' @param xinfo A list or matrix of cut points.
#' @return A list containing the model matrix, number of cut points, removed variables, and cut points information.
#' 
#' @export
#' 
#' @examples 
#' #' # Example usage
#' data(iris)
#' X <- iris[, -5]
#' numcut <- 3
#' usequants <- TRUE
#' rm.const <- TRUE
#' cont <- FALSE
#' xinfo <- NULL
#' 
#' #' result <- bartModelMatrix(X, numcut, usequants, rm.const, cont, xinfo)
#' #' print(result)
#' 
bartModelMatrix <- function(X,
                            numcut = 0L,
                            usequants = FALSE,
                            type = 7,
                            rm.const = FALSE,
                            cont = FALSE,
                            xinfo = NULL) {
    X.class <- class(X)[1]

    if (X.class == "factor") {
        X.class <- "data.frame"
        X <- data.frame(X = X)
    }

    grp <- NULL

    if (X.class == "data.frame") {
        p <- dim(X)[2]
        xnm <- names(X)
        for (i in 1:p) {
            if (is.factor(X[[i]])) {
                Xtemp <- class.ind(X[[i]])
                colnames(Xtemp) <- paste(xnm[i], 1:ncol(Xtemp), sep = "")
                X[[i]] <- Xtemp
                grp <- c(grp, rep(i, ncol(Xtemp)))
            } else {
                X[[i]] <- cbind(X[[i]])
                colnames(X[[i]]) <- xnm[i]
                grp <- c(grp, i)
            }
        }
        Xtemp <- cbind(X[[1]])
        if (p > 1) for (i in 2:p) Xtemp <- cbind(Xtemp, X[[i]])
        X <- Xtemp
    } else if (X.class == "numeric" | X.class == "integer") {
        X <- cbind(as.numeric(X))
        grp <- 1
    } else if (X.class == "NULL") {
        return(X)
    } else if (X.class != "matrix") {
        stop("Expecting either a factor, a vector, a matrix or a data.frame")
    }

    N <- nrow(X)
    p <- ncol(X)

    xinfo. <- matrix(nrow = p, ncol = numcut)
    nc <- numcut
    rm.vars <- c()

    if (N > 0 & p > 0 & (rm.const | numcut[1] > 0)) {
        for (j in 1:p) {
            X.class <- class(X[1, j])[1]

            if (X.class == "numeric" | X.class == "integer") {
                xs <- unique(sort(X[, j]))
                k <- length(xs)
                nc[j] <- numcut

                if (k %in% 0:1) {
                    rm.vars <- c(rm.vars, -j)
                    nc[j] <- 1
                    if (k == 0) xs <- NA
                } else if (cont) {
                    xs <- seq(xs[1], xs[k], length.out = numcut + 2)[-c(1, numcut + 2)]
                } else if (k < numcut) {
                    xs <- 0.5 * (xs[1:(k - 1)] + xs[2:k])
                    nc[j] <- k - 1
                } else if (usequants) {
                    xs <- quantile(X[, j],
                        type = type,
                        probs = (0:(numcut + 1)) / (numcut + 1)
                    )[-c(1, numcut + 2)]
                    names(xs) <- NULL
                } else {
                    xs <-
                        seq(xs[1], xs[k], length.out = numcut + 2)[-c(1, numcut + 2)]
                }
            } else {
                stop(paste0("Variables of type ", X.class, " are not supported"))
            }

            xinfo.[j, 1:nc[j]] <- xs
        }
    }

    X <- data.matrix(X)

    if (length(xinfo) > 0) {
        if (is.list(xinfo)) {
            for (j in 1:p) xinfo.[j, 1:length(xinfo[[j]])] <- xinfo[[j]]
        } else if (is.matrix(xinfo)) {
            xinfo. <- xinfo
        } else {
            stop("Only a list or a matrix can be provided for xinfo")
        }

        for (j in 1:p) nc[j] <- sum(!is.na(xinfo.[j, ]))
    }

    xinfo <- xinfo.

    if (rm.const & length(rm.vars) > 0) {
        X <- X[, rm.vars]
        nc <- nc[rm.vars]
        xinfo <- xinfo[rm.vars, ]
    } else if (length(rm.vars) == 0) rm.vars <- 1:p

    dimnames(xinfo) <- list(dimnames(X)[[2]], NULL)

    if (numcut == 0) {
        return(X)
    } else {
        return(list(
            X = X, numcut = as.integer(nc), rm.const = rm.vars,
            xinfo = xinfo, grp = grp
        ))
    }
}
