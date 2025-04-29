/*
 *  meBART: Bayesian Additive Regression Trees with Measurement Error
 *  Copyright (C) 2025 Kevin McCoy, Zachary Wooten, and Christine Peterson
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

#ifndef GUARD_mefuns_h
#define GUARD_mefuns_h

arma::mat rmvnorm(int n, arma::vec mu, arma::mat sigma);

arma::vec rmvnorm(arma::vec mu, arma::mat sigma);

double dmvnorm(arma::vec x, arma::vec mu, arma::mat sigma);

double rnorm(double mu, double sigma);

double dnorm(double x, double mu, double sigma);

#endif