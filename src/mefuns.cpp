/*
 * meBART: Bayesian Additive Regression Trees with Measurement Error
 * Copyright (C) 2025 Kevin McCoy, Zachary Wooten, and Christine Peterson
 * 
 * This package is a modification of the BART package originally created by
 * Robert McCulloch and Rodney Sparapani, who own the original copyright.
 * Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, a copy is available at
 * https://www.R-project.org/Licenses/GPL-2
 */

#include "common.h"
#include "mefuns.h"
#include <RcppArmadillo.h>
#include <random>

arn gen;

double rnorm(double mu, double sigma)
{
    return mu + sigma * gen.normal();
}

double dnorm(double x, double mu, double sigma)
{
    double tmp = (x - mu) / sigma;
    tmp = -0.5 * tmp * tmp;
    tmp = exp(tmp);
    tmp = tmp / (sigma * RTPI);
    return tmp;
}

arma::mat rmvnorm(int n, arma::vec mu, arma::mat sigma)
{
    int ncols = sigma.n_cols;
    arma::mat Y = arma::randn(n, ncols);
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

arma::vec rmvnorm(arma::vec mu, arma::mat sigma)
{
    int ncols = sigma.n_cols;
    arma::mat Y = arma::randn(1, ncols);
    arma::mat tmp = arma::repmat(mu, 1, 1).t() + Y * arma::chol(sigma);
    return tmp.row(0).t(); // return as a column vector
}

double dmvnorm(arma::vec x, arma::vec mu, arma::mat sigma)
{
    int dim = x.n_elem;
    double tmp = 1.0;

    tmp /= pow(RTPI, dim);
    tmp /= sqrt(arma::det(sigma));
    tmp *= exp(-0.5 * arma::as_scalar((x - mu).t() * arma::inv_sympd(sigma) * (x - mu)));

    return tmp;
}
