/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
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

#include <ctime>

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "rtnorm.h"
// #include "rtruncnorm.h"
#include "mefuns.h"

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)

RcppExport SEXP cpbart(
    SEXP _in,    // number of observations in training data
    SEXP _ip,    // dimension of x
    SEXP _inp,   // number of observations in test data
    SEXP _ix,    // x, train,  pxn (transposed so rows are contiguous in memory)
    SEXP _iy,    // y, train,  nx1
    SEXP _ixp,   // x, test, pxnp (transposed so rows are contiguous in memory)
    SEXP _im,    // number of trees
    SEXP _inc,   // number of cut points
    SEXP _ind,   // number of kept draws (except for thinnning ..)
    SEXP _iburn, // number of burn-in draws skipped
    SEXP _ipower,
    SEXP _ibase,
    SEXP _binaryOffset,
    SEXP _itau,
    //   SEXP _iM, // number of shards for Modified LISA MCMC
    SEXP _idart,
    SEXP _itheta,
    SEXP _iomega,
    SEXP _igrp,
    SEXP _ia,
    SEXP _ib,
    SEXP _irho,
    SEXP _iaug,
    SEXP _inkeeptrain,
    SEXP _inkeeptest,
    //   SEXP _inkeeptestme,
    SEXP _inkeeptreedraws,
    SEXP _inprintevery,
    //   SEXP _treesaslists,
    SEXP _Xinfo,           // cutpoints, now specified
    SEXP _proposal_sigma,  // standard deviation of proposal distribution for new MH step
    SEXP _meas_error_sigma, // standard deviation of measurement error
    SEXP _x_mu,
    SEXP _x_sigma // TODO: Add verbose parameter to control output
    )
{

    //--------------------------------------------------
    // process args
    size_t n = Rcpp::as<int>(_in);
    size_t p = Rcpp::as<int>(_ip);
    size_t np = Rcpp::as<int>(_inp);
    Rcpp::NumericVector xv(_ix);
    
    arma::mat x_matrix = Rcpp::as<arma::mat>(_ix); // me

    double *ix = &xv[0];
    Rcpp::IntegerVector yv(_iy); // binary
    int *iy = &yv[0];
    // Rcpp::NumericVector  xpv(_ixp);
    // double *ixp = &xpv[0];
    Rcpp::NumericMatrix xpv(_ixp);
    double *ixp = nullptr;
    if (np > 0)
        ixp = &xpv[0];
    size_t m = Rcpp::as<int>(_im);
    // size_t nc = Rcpp::as<int>(_inc);
    Rcpp::IntegerVector _nc(_inc);
    int *numcut = &_nc[0];
    size_t nd = Rcpp::as<int>(_ind);
    size_t burn = Rcpp::as<int>(_iburn);
    double mybeta = Rcpp::as<double>(_ipower);
    double alpha = Rcpp::as<double>(_ibase);
    double binaryOffset = Rcpp::as<double>(_binaryOffset);
    double tau = Rcpp::as<double>(_itau);
    //   double rootM = sqrt(Rcpp::as<double>(_iM));
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
    // Rcpp::IntegerVector _grp(_igrp);
    // int *grp = &_grp[0]; // unused
    size_t nkeeptrain = Rcpp::as<int>(_inkeeptrain);
    size_t nkeeptest = Rcpp::as<int>(_inkeeptest);
    //   size_t nkeeptestme = Rcpp::as<int>(_inkeeptestme);
    size_t nkeeptreedraws = Rcpp::as<int>(_inkeeptreedraws);
    size_t printevery = Rcpp::as<int>(_inprintevery);
    //   int treesaslists = Rcpp::as<int>(_treesaslists);
    Rcpp::NumericMatrix Xinfo(_Xinfo);
    // Rcpp::List Xinfo(_Xinfo);
    //   Rcpp::IntegerMatrix varcount(nkeeptreedraws, p);

    // return data structures (using Rcpp)
    /*
       Rcpp::NumericVector trmean(n); //train
       Rcpp::NumericVector temean(np);
    */
    Rcpp::NumericMatrix trdraw(nkeeptrain, n);
    Rcpp::NumericMatrix tedraw(nkeeptest, np);
    //   Rcpp::List list_of_lists(nkeeptreedraws*treesaslists);
    Rcpp::NumericMatrix varprb(nkeeptreedraws, p);
    Rcpp::IntegerMatrix varcnt(nkeeptreedraws, p);


    // Convert R sigma objets to arma matrices
    arma::mat proposal_sigma = Rcpp::as<arma::mat>(_proposal_sigma);
    arma::mat meas_error_sigma = Rcpp::as<arma::mat>(_meas_error_sigma);
    arma::vec x_mu = Rcpp::as<arma::vec>(_x_mu);
    arma::mat x_sigma = Rcpp::as<arma::mat>(_x_sigma);

    // random number generation
    arn gen;

    bart bm(m);

    if (Xinfo.size() > 0)
    {
        xinfo _xi;
        _xi.resize(p);
        for (size_t i = 0; i < p; i++)
        {
            _xi[i].resize(numcut[i]);
            // Rcpp::IntegerVector cutpts(Xinfo[i]);
            for (size_t j = 0; j < (size_t)numcut[i]; j++)
                _xi[i][j] = Xinfo(i, j);
        }
        bm.setxinfo(_xi);
    }
#else

#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]

void cpbart(
    size_t n,    // number of observations in training data
    size_t p,    // dimension of x
    size_t np,   // number of observations in test data
    double *ix,  // x, train,  pxn (transposed so rows are contiguous in memory)
    int *iy,     // y, train,  nx1
    double *ixp, // x, test, pxnp (transposed so rows are contiguous in memory)
    size_t m,    // number of trees
    int *numcut, // size_t nc,		//number of cut points
    size_t nd,   // number of kept draws (except for thinnning ..)
    size_t burn, // number of burn-in draws skipped
    double mybeta,
    double alpha,
    double binaryOffset,
    double tau,
    bool dart,
    double theta,
    double omega,
    int *grp,
    double a,
    double b,
    double rho,
    bool aug,
    size_t nkeeptrain,
    size_t nkeeptest,
    // size_t nkeeptestme,
    size_t nkeeptreedraws,
    size_t printevery,
    //   int treesaslists,
    unsigned int n1, // additional parameters needed to call from C++
    unsigned int n2,
    /*
       double* trmean,
       double* temean,
    */
    double *_trdraw,
    double *_tedraw)
{

    // return data structures (using C++)
    std::vector<double *> trdraw(nkeeptrain);
    std::vector<double *> tedraw(nkeeptest);

    for (size_t i = 0; i < nkeeptrain; ++i)
        trdraw[i] = &_trdraw[i * n];
    for (size_t i = 0; i < nkeeptest; ++i)
        tedraw[i] = &_tedraw[i * np];

    std::vector<std::vector<size_t>> varcnt;
    std::vector<std::vector<double>> varprb;

    // random number generation
    arn gen(n1, n2);

    bart bm(m);
#endif

    /*
       for(size_t i=0; i<n; ++i) trmean[i]=0.0;
       for(size_t i=0; i<np; ++i) temean[i]=0.0;
    */

    double *iz = new double[n];

    std::stringstream treess; // string stream to write trees to
    treess.precision(10);
    treess << nkeeptreedraws << " " << m << " " << p << endl;
    // dart iterations
    std::vector<double> ivarprb(p, 0.);
    std::vector<size_t> ivarcnt(p, 0);

    printf("*****Into main of pbart\n");

    size_t skiptr, skipte, skiptreedraws;
    // size_t skiptr,skipte,skipteme,skiptreedraws;
    if (nkeeptrain)
    {
        skiptr = nd / nkeeptrain;
    }
    else
        skiptr = nd + 1;
    if (nkeeptest)
    {
        skipte = nd / nkeeptest;
    }
    else
        skipte = nd + 1;
    /*
       if(nkeeptestme) {skipteme=nd/nkeeptestme;}
       else skipteme=nd+1;
    */
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
    printf("y1,yn: %d, %d\n", iy[0], iy[n - 1]);
    printf("x1,x[n*p]: %lf, %lf\n", ix[0], ix[n * p - 1]);
    if (np)
        printf("xp1,xp[np*p]: %lf, %lf\n", ixp[0], ixp[np * p - 1]);
    printf("*****Number of Trees: %zu\n", m);
    printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p - 1]);
    printf("*****burn and ndpost: %zu, %zu\n", burn, nd);
    printf("*****Prior:mybeta,alpha,tau: %lf,%lf,%lf\n",
           mybeta, alpha, tau);
    printf("*****binaryOffset: %lf\n", binaryOffset);
    cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: "
         << dart << ',' << theta << ',' << omega << ',' << a << ','
         << b << ',' << rho << ',' << aug << endl;
    printf("*****nkeeptrain,nkeeptest,nkeeptreedraws: %zu,%zu,%zu\n",
           nkeeptrain, nkeeptest, nkeeptreedraws);
    /*
       printf("*****nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws: %zu,%zu,%zu,%zu\n",
                   nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws);
    */
    printf("*****printevery: %zu\n", printevery);
    // printf("*****skiptr,skipte,skipteme,skiptreedraws: %zu,%zu,%zu,%zu\n",skiptr,skipte,skipteme,skiptreedraws);
    printf("*****skiptr,skipte,skiptreedraws: %zu,%zu,%zu\n", skiptr, skipte, skiptreedraws);

    //--------------------------------------------------
    // bart bm(m);
    bm.setprior(alpha, mybeta, tau);
    bm.setdata(p, n, ix, iz, numcut);
    bm.setdart(a, b, rho, aug, dart, theta, omega);
    //--------------------------------------------------
    // init
    for (size_t k = 0; k < n; k++)
    {
        if (iy[k] == 0)
            iz[k] = -rtnorm(0., binaryOffset, 1., gen);
        else
            iz[k] = rtnorm(0., -binaryOffset, 1., gen);
        /*
            if(iy[k]==0) iz[k]= -r_lefttruncnorm(0., binaryOffset, 1., gen);
            else iz[k]=r_lefttruncnorm(0., -binaryOffset, 1., gen);
        */
    }

    //--------------------------------------------------

    //--------------------------------------------------
    // temporary storage
    // out of sample fit
    double *fhattest = 0;
    if (np)
    {
        fhattest = new double[np];
    }

    //--------------------------------------------------
    // mcmc
    printf("\nMCMC\n");
    // size_t index;
    size_t trcnt = 0; // count kept train draws
    size_t tecnt = 0; // count kept test draws
    //   size_t temecnt=0; //count test draws into posterior mean
    //   size_t treedrawscnt=0; //count kept bart draws
    bool keeptest, /*keeptestme,*/ keeptreedraw;

    time_t tp;
    int time1 = time(&tp);
    xinfo &xi = bm.getxinfo();
    size_t total = nd + burn;


    // Define x_draws object to store all draws of x
    arma::cube x_draws_(p, n, total + 1); // p x n x (total + 1) cube to store all draws of x
    x_draws_.slice(0) = x_matrix; // Initialize first entry of x_draws with the observed x values

    // Create storage for acceptances of each x draw
    Rcpp::NumericMatrix acceptances(total, n);



    for (size_t i = 0; i < total; i++)
    {
        if (i % printevery == 0)
            printf("done %zu (out of %zu)\n", i, total);
        if (i == (burn / 2) && dart)
            bm.startdart();
        // draw bart
        bm.draw(1., gen);
        // bm.draw(rootM, gen);

        for (size_t k = 0; k < n; k++)
        {
            if (iy[k] == 0)
                iz[k] = -rtnorm(-bm.f(k), binaryOffset, 1., gen);
            else
                iz[k] = rtnorm(bm.f(k), -binaryOffset, 1., gen);
            /*
                if(iy[k]==0) iz[k]= -r_lefttruncnorm(-bm.f(k), binaryOffset, 1., gen);
                else iz[k]=r_lefttruncnorm(bm.f(k), -binaryOffset, 1., gen);
            */
        }

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


            double y_true = iz[k];
            double y_pred = bm.f(k);


            // TODO: Check that this is right, probably isnt
            bart bm_prime;
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
            bm_prime.setdata(p, n, ix, iz, numcut);
            double y_pred_prime = bm_prime.f(k);

            // Calculate MH ratio
            double alpha = 0.0;

            // Old values
            alpha -= log(dnorm(y_true, y_pred, 1.)); // y likelihood
            alpha -= log(dmvnorm(x_meas, x_true, meas_error_sigma)); // x likelihood
            alpha -= log(dmvnorm(x_true, x_mu, x_sigma));   // x prior

            // Proposed values
            alpha += log(dnorm(y_true, y_pred_prime, 1.));   // y likelihood
            alpha += log(dmvnorm(x_meas, x_true_prime, meas_error_sigma)); // x likelihood
            alpha += log(dmvnorm(x_true_prime, x_mu, x_sigma));   // x prior


            // Calculate Metropolis-Hastings acceptance ratio
            alpha = exp(alpha); // Convert back from log scale
            double acceptance_ratio = std::min<double>(1, alpha);
            bool accept = gen.uniform() < acceptance_ratio;

            // Update draw of X
            if (accept)
            {
                x_draws_.slice(i + 1).col(k) = x_true_prime; // Update the draw of x_true

                // Give original BART model new x data
                bm.resetdata(p, n, ix, iz);
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

            // for(size_t k=0;k<n;k++) trmean[k]+=bm.f(k);
            if (nkeeptrain && (((i - burn + 1) % skiptr) == 0))
            {
                // index = trcnt*n;;
                // for(size_t k=0;k<n;k++) trdraw[index+k]=bm.f(k);
                for (size_t k = 0; k < n; k++)
                    TRDRAW(trcnt, k) = bm.f(k);
                trcnt += 1;
            }
            keeptest = nkeeptest && (((i - burn + 1) % skipte) == 0) && np;
            // keeptestme = nkeeptestme && (((i-burn+1) % skipteme) ==0) && np;
            // if(keeptest || keeptestme) bm.predict(p,np,ixp,fhattest);
            if (keeptest)
            {
                bm.predict(p, np, ixp, fhattest);
                // index=tecnt*np;
                // for(size_t k=0;k<np;k++) tedraw[index+k]=fhattest[k];
                for (size_t k = 0; k < np; k++)
                    TEDRAW(tecnt, k) = fhattest[k];
                tecnt += 1;
            }
            /*
                     if(keeptestme) {
                        for(size_t k=0;k<np;k++) temean[k]+=fhattest[k];
                        temecnt+=1;
                     }
            */
            keeptreedraw = nkeeptreedraws && (((i - burn + 1) % skiptreedraws) == 0);
            if (keeptreedraw)
            {
                //	   #ifndef NoRcpp
                //	   Rcpp::List lists(m*treesaslists);
                //	   #endif

                for (size_t j = 0; j < m; j++)
                {
                    treess << bm.gettree(j);
                    /*
                              #ifndef NoRcpp
                              varcount.row(treedrawscnt)=varcount.row(treedrawscnt)+bm.gettree(j).tree2count(p);
                              if(treesaslists) lists(j)=bm.gettree(j).tree2list(xi, 0., 1.);
                              #endif
                    */
                }
#ifndef NoRcpp
                //	    if(treesaslists) list_of_lists(treedrawscnt)=lists;
                ivarcnt = bm.getnv();
                ivarprb = bm.getpv();
                size_t k = (i - burn) / skiptreedraws;
                for (size_t j = 0; j < p; j++)
                {
                    varcnt(k, j) = ivarcnt[j];
                    // varcnt(i-burn,j)=ivarcnt[j];
                    varprb(k, j) = ivarprb[j];
                    // varprb(i-burn,j)=ivarprb[j];
                }
#else
                varcnt.push_back(bm.getnv());
                varprb.push_back(bm.getpv());
#endif

                // treedrawscnt +=1;
            }
        }
    }
    int time2 = time(&tp);
    printf("time: %ds\n", time2 - time1);
    /*
       for(size_t k=0;k<n;k++) trmean[k]/=nd;
       for(size_t k=0;k<np;k++) temean[k]/=temecnt;
    */
    printf("check counts\n");
    printf("trcnt,tecnt: %zu,%zu\n", trcnt, tecnt);
    // printf("trcnt,tecnt,temecnt,treedrawscnt: %zu,%zu,%zu,%zu\n",trcnt,tecnt,temecnt,treedrawscnt);
    //--------------------------------------------------
    // PutRNGstate();

    if (fhattest)
        delete[] fhattest;
    delete[] iz;

#ifndef NoRcpp
    //--------------------------------------------------
    // return
    Rcpp::List ret;
    // ret["yhat.train.mean"]=trmean;
    ret["yhat.train"] = trdraw;
    // ret["yhat.test.mean"]=temean;
    ret["yhat.test"] = tedraw;
    // ret["varcount"]=varcount;
    ret["varcount"] = varcnt;
    ret["varprob"] = varprb;

    // for(size_t i=0;i<m;i++) {
    //   bm.gettree(i).pr();
    //}

    Rcpp::List xiret(xi.size());
    for (size_t i = 0; i < xi.size(); i++)
    {
        Rcpp::NumericVector vtemp(xi[i].size());
        std::copy(xi[i].begin(), xi[i].end(), vtemp.begin());
        xiret[i] = Rcpp::NumericVector(vtemp);
    }

    Rcpp::List treesL;
    // treesL["nkeeptreedraws"] = Rcpp::wrap<int>(nkeeptreedraws); //in trees
    // treesL["ntree"] = Rcpp::wrap<int>(m); //in trees
    // treesL["numx"] = Rcpp::wrap<int>(p); //in cutpoints
    treesL["cutpoints"] = xiret;
    treesL["trees"] = Rcpp::CharacterVector(treess.str());
    // if(treesaslists) treesL["lists"]=list_of_lists;
    ret["treedraws"] = treesL;

    ret["x_draws"] = x_draws_; // TODO: Don't return burn-in values
    ret["acceptances"] = acceptances; // TODO: Don't return burn-in values

    return ret;
#else

#endif
}