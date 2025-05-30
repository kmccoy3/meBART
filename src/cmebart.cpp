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

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "heterbart.h"
#include "mefuns.h"

#define TRDRAW(a, b) trdraw(a, b) // pre-processor directive to replace former with latter
#define TEDRAW(a, b) tedraw(a, b) // same as above

RcppExport SEXP cmebart(
    SEXP _in,              // number of observations in training data
    SEXP _ip,              // dimension of x
    SEXP _inp,             // number of observations in test data
    SEXP _ix,              // x, train,  pxn (transposed so rows are contiguous in memory)
    SEXP _iy,              // y, train,  nx1
    SEXP _ixp,             // x, test, pxnp (transposed so rows are contiguous in memory)
    SEXP _im,              // number of trees
    SEXP _inc,             // number of cut points
    SEXP _ind,             // number of kept draws (except for thinnning ..)
    SEXP _iburn,           // number of burn-in draws skipped
    SEXP _ipower,          // tree depth prior, beta in CGM
    SEXP _ibase,           // tree depth prior, alpha in CGM
    SEXP _itau,            // sigma from mu prior, now scaled by number of trees
    SEXP _inu,             // DOF in sigma prior
    SEXP _ilambda,         // scale parameter used in sigma prior
    SEXP _isigest,         // sigma prior, estimated from the data
    SEXP _iw,              // vector of weights to multiply the standard deviation, defaults to one
    SEXP _idart,           // use Dirichlet prior for variable selection
    SEXP _itheta,          // TODO ???
    SEXP _iomega,          // TODO ???
    SEXP _igrp,            // counts out which variables are unique / which ones are factors
    SEXP _ia,              // Beta parameters, used in DART only
    SEXP _ib,              // Beta parameters, used in DART only
    SEXP _irho,            // Sparsity parameter, used in DART only
    SEXP _iaug,            // TODO ???
    SEXP _inkeeptrain,     // number of training data posterior draws to keep
    SEXP _inkeeptest,      // number of test data posterior draws to keep
    SEXP _inkeeptestme,    // number of MCMC iters to be returned for test mean
    SEXP _inkeeptreedraws, // number of MCMC iters to be returned for tree draws
    SEXP _inprintevery,    // print progress every printevery iterations
                           //   SEXP _treesaslists,
    SEXP _Xinfo,           // cutpoints, now specified
    SEXP _proposal_sigma,  // standard deviation of proposal distribution for new MH step
    SEXP _meas_error_sigma, // standard deviation of measurement error
    SEXP _x_mu,
    SEXP _x_sigma // TODO: Add verbose parameter to control output
)
{
    // printf("FIRST HEADER RAN\n");
    //--------------------------------------------------
    // process args
    size_t n = Rcpp::as<int>(_in);   // number of training data points
    size_t p = Rcpp::as<int>(_ip);   // number of predictors
    size_t np = Rcpp::as<int>(_inp); // number of test data points


    Rcpp::NumericVector xv(_ix);
    arma::mat x_matrix = Rcpp::as<arma::mat>(_ix);


    double *ix = &xv[0];
    Rcpp::NumericVector yv(_iy);
    double *iy = &yv[0];

    Rcpp::NumericMatrix xpv(_ixp); // test data
    double *ixp = nullptr;
    if (np > 0)
        ixp = &xpv[0];
    size_t m = Rcpp::as<int>(_im); // number of trees
    Rcpp::IntegerVector _nc(_inc); // number of cut points
    int *numcut = &_nc[0];
    size_t nd = Rcpp::as<int>(_ind);
    size_t burn = Rcpp::as<int>(_iburn);
    double mybeta = Rcpp::as<double>(_ipower);
    double alpha = Rcpp::as<double>(_ibase);
    double tau = Rcpp::as<double>(_itau);
    double nu = Rcpp::as<double>(_inu);
    double lambda = Rcpp::as<double>(_ilambda);
    double sigma = Rcpp::as<double>(_isigest);
    Rcpp::NumericVector wv(_iw);
    double *iw = &wv[0];
    bool dart;
    if (Rcpp::as<int>(_idart) == 1)
        dart = true;
    else
        dart = false;
    double a = Rcpp::as<double>(_ia);
    double b = Rcpp::as<double>(_ib);
    double rho = Rcpp::as<double>(_irho);
    bool aug;
    if (Rcpp::as<int>(_iaug) == 1)
        aug = true;
    else
        aug = false;
    double theta = Rcpp::as<double>(_itheta);
    double omega = Rcpp::as<double>(_iomega);
    // Rcpp::IntegerVector _grp(_igrp); // I removed this, was causing warning because unused
    // int *grp = &_grp[0];
    size_t nkeeptrain = Rcpp::as<int>(_inkeeptrain);
    size_t nkeeptest = Rcpp::as<int>(_inkeeptest);
    size_t nkeeptestme = Rcpp::as<int>(_inkeeptestme);
    size_t nkeeptreedraws = Rcpp::as<int>(_inkeeptreedraws);
    size_t printevery = Rcpp::as<int>(_inprintevery);
    Rcpp::NumericMatrix Xinfo(_Xinfo);

    // return data structures (using Rcpp)
    Rcpp::NumericVector trmean(n); // train
    Rcpp::NumericVector temean(np);
    Rcpp::NumericVector sdraw(nd + burn);
    Rcpp::NumericMatrix trdraw(nkeeptrain, n);
    Rcpp::NumericMatrix tedraw(nkeeptest, np);
    Rcpp::NumericMatrix varprb(nkeeptreedraws, p);
    Rcpp::IntegerMatrix varcnt(nkeeptreedraws, p);

    // Convert R sigma objets to arma matrices
    arma::mat proposal_sigma = Rcpp::as<arma::mat>(_proposal_sigma);
    arma::mat meas_error_sigma = Rcpp::as<arma::mat>(_meas_error_sigma);
    arma::vec x_mu = Rcpp::as<arma::vec>(_x_mu);
    arma::mat x_sigma = Rcpp::as<arma::mat>(_x_sigma);

    // random number generation
    arn gen; // TODO

    heterbart bm(m); // Initializes a heterbart object with m trees

    //-----------------------------------------------------------
    // This section really just transforms the x cut points, makes it part of the bm object
    if (Xinfo.size() > 0)
    {
        // typedef std::vector<vec_d> xinfo; //vector of vectors, will be split rules
        xinfo _xi;
        _xi.resize(p); // Makes outside vector p long
        for (size_t i = 0; i < p; i++)
        {
            _xi[i].resize(numcut[i]); // resizes each internal vector to correct number of cut points
            for (size_t j = 0; j < (size_t)numcut[i]; j++)
                _xi[i][j] = Xinfo(i, j); // Sets values to the cut points
        }
        bm.setxinfo(_xi);
    }
    //-----------------------------------------------------------
    // initializes the vectors to zero

    for (size_t i = 0; i < n; i++)
        trmean[i] = 0.0;
    for (size_t i = 0; i < np; i++)
        temean[i] = 0.0;

    printf("*****Into main of wbart --- KEVINS CODE :)\n");
    //-----------------------------------------------------------
    // Calculates the skip values for the MCMC iterations
    size_t skiptr, skipte, skipteme, skiptreedraws;
    if (nkeeptrain)
    {
        skiptr = nd / nkeeptrain; // How many iters to skip during training data
    }
    else
        skiptr = nd + 1;
    if (nkeeptest)
    {
        skipte = nd / nkeeptest;
    }
    else
        skipte = nd + 1; // How many iters to skip during test data
    if (nkeeptestme)
    {
        skipteme = nd / nkeeptestme; // How many iters to skip during test data for posterior mean
    }
    else
        skipteme = nd + 1;
    if (nkeeptreedraws)
    {
        skiptreedraws = nd / nkeeptreedraws;
    }
    else
        skiptreedraws = nd + 1;

    //--------------------------------------------------
    // print args
    printf("*****Data:\n");
    printf("data:n,p,np: %zu, %zu, %zu\n", n, p, np);
    printf("y1,yn: %lf, %lf\n", iy[0], iy[n - 1]);
    printf("x1,x[n*p]: %lf, %lf\n", ix[0], ix[n * p - 1]);
    if (np)
        printf("xp1,xp[np*p]: %lf, %lf\n", ixp[0], ixp[np * p - 1]);
    printf("*****Number of Trees: %zu\n", m);
    printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p - 1]);
    printf("*****burn and ndpost: %zu, %zu\n", burn, nd);
    printf("*****Prior:beta,alpha,tau,nu,lambda: %lf,%lf,%lf,%lf,%lf\n",
           mybeta, alpha, tau, nu, lambda);
    printf("*****sigma: %lf\n", sigma);
    printf("*****w (weights): %lf ... %lf\n", iw[0], iw[n - 1]);
    cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: "
         << dart << ',' << theta << ',' << omega << ',' << a << ','
         << b << ',' << rho << ',' << aug << endl;
    printf("*****nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws: %zu,%zu,%zu,%zu\n",
           nkeeptrain, nkeeptest, nkeeptestme, nkeeptreedraws);
    printf("*****printevery: %zu\n", printevery);
    printf("*****skiptr,skipte,skipteme,skiptreedraws: %zu,%zu,%zu,%zu\n", skiptr, skipte, skipteme, skiptreedraws);

    //--------------------------------------------------

    bm.setprior(alpha, mybeta, tau);
    bm.setdata(p, n, ix, iy, numcut); // TODO: makes data part of the bm object
    bm.setdart(a, b, rho, aug, dart, theta, omega);

    //--------------------------------------------------
    // sigma = sigest, adjusts prior of mu
    double *svec = new double[n];
    for (size_t i = 0; i < n; i++)
        svec[i] = iw[i] * sigma; // incorporate weights into sigma hat estimate

    //--------------------------------------------------

    std::stringstream treess; // string stream to write trees to
    treess.precision(10);
    treess << nkeeptreedraws << " " << m << " " << p << endl;
    // dart iterations
    std::vector<double> ivarprb(p, 0.);
    std::vector<size_t> ivarcnt(p, 0);

    //--------------------------------------------------
    // temporary storage
    // out of sample fit
    double *fhattest = 0; // posterior mean for prediction
    if (np)
    {
        fhattest = new double[np];
    }
    double restemp = 0.0, rss = 0.0;

    //--------------------------------------------------
    // mcmc
    printf("\nMCMC  -------- KEVINS CODE :)\n");

    size_t trcnt = 0;        // count kept train draws
    size_t tecnt = 0;        // count kept test draws
    size_t temecnt = 0;      // count test draws into posterior mean
    size_t treedrawscnt = 0; // count kept bart draws
    bool keeptest, keeptestme, keeptreedraw;

    time_t tp;
    int time1 = time(&tp);     // start time
    xinfo &xi = bm.getxinfo(); // get cutpoints back
    size_t total = nd + burn;  // total number of MCMC iterations

    // Define x_draws object to store all draws of x
    arma::cube x_draws_(p, n, total + 1); // p x n x (total + 1) cube to store all draws of x
    x_draws_.slice(0) = x_matrix; // Initialize first entry of x_draws with the observed x values

    // Create storage for acceptances of each x draw
    Rcpp::NumericMatrix acceptances(total, n);

    for (size_t i = 0; i < total; i++) // main MCMC loop
    {
        if (i % printevery == 0)
            printf("done %zu (out of %zu)\n", i, total);
        if (i == (burn / 2) && dart)
            bm.startdart();
        // draw bart
        bm.draw(svec, gen);

        // Calculate the residual sum of squares
        rss = 0.0;
        for (size_t k = 0; k < n; k++)
        {
            restemp = (iy[k] - bm.f(k)) / (iw[k]); // residual = (y - f(x)) / w
            rss += restemp * restemp;
        }
        // Calculate sigma / draw sigma
        sigma = sqrt((nu * lambda + rss) / gen.chi_square(n + nu));
        // Recalculate svec
        for (size_t k = 0; k < n; k++)
            svec[k] = iw[k] * sigma;
        // Save sigma to sdraw
        sdraw[i] = sigma;

        // =========================================================================================
        // Measurement Error Step in Gibbs Sampler

        Rcpp::NumericVector last_xv(n*p);
        for (size_t j = 0; j < p ; j++){
            for (size_t l = 0; l < n; l++){
                last_xv[j + l*p] = x_draws_.slice(i).col(l)[j]; // Copy the current x_draws_ into last_x
            }
        }


        // Loop through each observation
        for (size_t k = 0; k < n; k++)
        {

            // Get x values of interest
            arma::vec x_meas = x_matrix.col(k); // observed value of x
            arma::vec x_true = x_draws_.slice(i).col(k); // old value of x_true
            arma::vec x_true_prime = rmvnorm(x_true, proposal_sigma); // // TODO: Fix hardcoding of 0.1


            //     // Hyperparameters
            //     // double mu_x = 0.5;     // Prior mean
            //     // double sigma_x = 0.25; // Prior standard deviation
            //     double sigma_e = 0.1; // Measurement error standard deviation

            //     // Initialize alpha

            double y_true = yv[k];
            double y_pred = bm.f(k);
            //     double y_pred_prime; // = y_pred; // FIXME: we need to calculate y_pred for the new x_true_prime

            // TODO: Check that this is right, probably isnt
            heterbart bm_prime;
            bm_prime = bm;


            Rcpp::NumericVector xv_prime(Rcpp::clone(last_xv));
            for (size_t j=0; j<p; j++){
                xv_prime[j + k*p] = x_true_prime[j];
            }

            // if (i < 5){
            //     Rcpp::Rcout << "x_meas: " << xv << std::endl;
            //     Rcpp::Rcout << "last_xv: " << last_xv << std::endl;
            // }

            double *ix = &xv_prime[0];
            // double *ix = &last_xv[0];
            bm_prime.setdata(p, n, ix, iy, numcut);
            double y_pred_prime = bm_prime.f(k);

            // Calculate MH ratio
            double alpha = 0.0;

            // Old values
            alpha -= log(dnorm(y_true, y_pred, sigma)); // y likelihood
            alpha -= log(dmvnorm(x_meas, x_true, meas_error_sigma)); // x likelihood
            alpha -= log(dmvnorm(x_true, x_mu, x_sigma));   // x prior

            // Proposed values
            alpha += log(dnorm(y_true, y_pred_prime, sigma));   // y likelihood
            alpha += log(dmvnorm(x_meas, x_true_prime, meas_error_sigma)); // x likelihood
            alpha += log(dmvnorm(x_true_prime, x_mu, x_sigma));   // x prior


            // Calculate Metropolis-Hastings acceptance ratio
            alpha = exp(alpha); // Convert back from log scale
            double acceptance_ratio = std::min<double>(1, alpha);
            bool accept = gen.uniform() < acceptance_ratio;

            // Update draw of X
            if (accept && i > burn)
            {
                x_draws_.slice(i + 1).col(k) = x_true_prime; // Update the draw of x_true

                // Give original BART model new x data
                bm.resetdata(p, n, ix, iy);
                acceptances(i, k) = 1;
            }
            else
            {
                x_draws_.slice(i + 1).col(k) = x_true; // Keep the old value of x_true
                acceptances(i, k) = 0;
            }
        }

        // =========================================================================================

        if (i >= burn)
        {
            // add to the mean of the training data
            for (size_t k = 0; k < n; k++)
                trmean[k] += bm.f(k);

            // if ...
            if (nkeeptrain && (((i - burn + 1) % skiptr) == 0))
            {

                for (size_t k = 0; k < n; k++)
                    TRDRAW(trcnt, k) = bm.f(k);
                trcnt += 1;
            }
            keeptest = nkeeptest && (((i - burn + 1) % skipte) == 0) && np;       // if we are keeping test data
            keeptestme = nkeeptestme && (((i - burn + 1) % skipteme) == 0) && np; // if we are keeping test data for posterior mean
            if (keeptest || keeptestme)
                bm.predict(p, np, ixp, fhattest);
            if (keeptest)
            {

                for (size_t k = 0; k < np; k++)
                    TEDRAW(tecnt, k) = fhattest[k];
                tecnt += 1;
            }
            if (keeptestme)
            {
                for (size_t k = 0; k < np; k++)
                    temean[k] += fhattest[k];
                temecnt += 1;
            }
            keeptreedraw = nkeeptreedraws && (((i - burn + 1) % skiptreedraws) == 0);
            if (keeptreedraw)
            {

                for (size_t j = 0; j < m; j++)
                {
                    treess << bm.gettree(j); // write tree to string stream
                }

                ivarcnt = bm.getnv();
                ivarprb = bm.getpv();
                size_t k = (i - burn) / skiptreedraws;
                for (size_t j = 0; j < p; j++)
                {
                    varcnt(k, j) = ivarcnt[j];

                    varprb(k, j) = ivarprb[j];
                }

                treedrawscnt += 1;
            }
        }
    } // end of MCMC loop
    //--------------------------------------------------

    int time2 = time(&tp);
    printf("time: %ds\n", time2 - time1);
    for (size_t k = 0; k < n; k++)
        trmean[k] /= nd;
    for (size_t k = 0; k < np; k++)
        temean[k] /= temecnt;
    printf("check counts\n");
    printf("trcnt,tecnt,temecnt,treedrawscnt: %zu,%zu,%zu,%zu\n", trcnt, tecnt, temecnt, treedrawscnt);
    //--------------------------------------------------

    if (fhattest)
        delete[] fhattest;
    if (svec)
        delete[] svec;

    //--------------------------------------------------
    // return

    Rcpp::List ret;
    ret["sigma"] = sdraw;
    ret["yhat.train.mean"] = trmean;
    ret["yhat.train"] = trdraw;
    ret["yhat.test.mean"] = temean;
    ret["yhat.test"] = tedraw;
    ret["varcount"] = varcnt;
    ret["varprob"] = varprb;

    Rcpp::List xiret(xi.size());
    for (size_t i = 0; i < xi.size(); i++)
    {
        Rcpp::NumericVector vtemp(xi[i].size());
        std::copy(xi[i].begin(), xi[i].end(), vtemp.begin());
        xiret[i] = Rcpp::NumericVector(vtemp);
    }

    Rcpp::List treesL;
    treesL["cutpoints"] = xiret;
    treesL["trees"] = Rcpp::CharacterVector(treess.str());
    ret["treedraws"] = treesL;

    ret["x_draws"] = x_draws_; // TODO: Don't return burn-in values
    ret["acceptances"] = acceptances; // TODO: Don't return burn-in values

    return ret;
}
