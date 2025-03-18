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

#include "bart.h"

//--------------------------------------------------
// constructor
bart::bart() : m(200), t(m), pi(), p(0), n(0), x(0), y(0), xi(), allfit(0), r(0), ftemp(0), di(), dartOn(false) {}
bart::bart(size_t im) : m(im), t(m), pi(), p(0), n(0), x(0), y(0), xi(), allfit(0), r(0), ftemp(0), di(), dartOn(false) {}
bart::bart(const bart &ib) : m(ib.m), t(m), pi(ib.pi), p(0), n(0), x(0), y(0), xi(), allfit(0), r(0), ftemp(0), di(), dartOn(false)
{
    this->t = ib.t;
}
bart::~bart()
{
    if (allfit)
        delete[] allfit;
    if (r)
        delete[] r;
    if (ftemp)
        delete[] ftemp;
}

//--------------------------------------------------
// operators
bart &bart::operator=(const bart &rhs)
{
    if (&rhs != this)
    {

        this->t = rhs.t;
        this->m = t.size();

        this->pi = rhs.pi;

        p = 0;
        n = 0;
        x = 0;
        y = 0;
        xi.clear();

        if (allfit)
        {
            delete[] allfit;
            allfit = 0;
        }
        if (r)
        {
            delete[] r;
            r = 0;
        }
        if (ftemp)
        {
            delete[] ftemp;
            ftemp = 0;
        }
    }
    return *this;
}
//--------------------------------------------------
// get,set
void bart::setm(size_t m)
{
    t.resize(m); // resize vector of trees
    this->m = t.size(); // set m

    if (allfit && (xi.size() == p)) // if allfit exists and xi is the same size as p
        predict(p, n, x, allfit);
}

//--------------------------------------------------
void bart::setxinfo(xinfo &_xi)
{
    size_t p = _xi.size();
    xi.resize(p);
    for (size_t i = 0; i < p; i++)
    {
        size_t nc = _xi[i].size();
        xi[i].resize(nc);
        for (size_t j = 0; j < nc; j++)
            xi[i][j] = _xi[i][j];
    }
}
//--------------------------------------------------
// this version is when numcut is the same for all variables, so DOESNT RUN
void bart::setdata(size_t p, size_t n, double *x, double *y, size_t numcut)
// dimension p, size n, pointer to x VECTOR, pointer to y vector, number of cutpoints
{
    // printf("In bart::setdata!!!!!!!!!\n");
    int *nc = new int[p]; // number of cutpoints for each variable, p-length vector
    for (size_t i = 0; i < p; ++i)
        nc[i] = numcut; // loop over each variable, set number of cutpoints to numcut
    this->setdata(p, n, x, y, nc);
    delete[] nc;
}

void bart::setdata(size_t p, size_t n, double *x, double *y, int *nc)
{
    // printf("In bart::setdata v2!!!!!!!!!\n");
    this->p = p;
    this->n = n;
    this->x = x;

    // printf("Address of x is: %p\n", x);
    // printf("The value of x is: %f\n", *x);

    this->y = y; // sets inputs as member variables

    if (xi.size() == 0) // if xi (vector of vectors, cutpoints) is empty
        makexinfo(p, n, &x[0], xi, nc);

    if (allfit)
        delete[] allfit;
    allfit = new double[n];
    predict(p, n, x, allfit);

    if (r)
        delete[] r;
    r = new double[n];

    if (ftemp)
        delete[] ftemp;
    ftemp = new double[n];

    di.n = n;
    di.p = p;
    di.x = &x[0];
    di.y = r;
    for (size_t j = 0; j < p; j++)
    {
        nv.push_back(0);
        pv.push_back(1 / (double)p);
    }
}


void bart::resetdata(size_t p, size_t n, double *x, double *y)
{
    // printf("In bart::resetdata!!!\n");
    // this->p = p;
    // this->n = n;

    // printf("Reassigning x\n");
    this->x = x;

    // printf("Address of x is: %p\n", x);
    // printf("The value of x is: %f\n", *x);

    // this->y = y; // sets inputs as member variables

    // if (xi.size() == 0) // if xi (vector of vectors, cutpoints) is empty
    //     makexinfo(p, n, &x[0], xi, nc);

    // printf("Resetting allfit\n");
    if (allfit)
        delete[] allfit;
    allfit = new double[n];
    predict(p, n, x, allfit);


    // printf("Resetting r\n");
    // if (r)
    //     delete[] r;
    // r = new double[n];


    // printf("Resetting ftemp\n");
    // if (ftemp)
    //     delete[] ftemp;
    // ftemp = new double[n];


    // printf("Resetting di\n");
    // di.n = n;
    // di.p = p;
    di.x = &x[0];
    // di.y = r;
    // for (size_t j = 0; j < p; j++)
    // {
    //     nv.push_back(0);
    //     pv.push_back(1 / (double)p);
    // }
}





//--------------------------------------------------
void bart::predict(size_t p, size_t n, double *x, double *fp)
// uses: m,t,xi
{
    double *fptemp = new double[n];

    for (size_t j = 0; j < n; j++) // loop over observations
        fp[j] = 0.0;               // this is the allfit object, n-vector or now floats

    for (size_t j = 0; j < m; j++) // loop over trees
    {
        fit(t[j], xi, p, n, x, fptemp); // fit the tree to the data
        for (size_t k = 0; k < n; k++)  // loop over observations
            fp[k] += fptemp[k];         // add the tree's prediction to the overall prediction
    }

    delete[] fptemp;
}
//--------------------------------------------------
void bart::draw(double sigma, rn &gen)
{
    // Loop over trees
    for (size_t j = 0; j < m; j++)
    {
        fit(t[j], xi, p, n, x, ftemp); // fit the tree to the data
        // Loop over observations
        for (size_t k = 0; k < n; k++)
        {
            allfit[k] = allfit[k] - ftemp[k]; // takes out the current tree's prediction
            r[k] = y[k] - allfit[k]; // recalculates residuals
        }
        bd(t[j], xi, di, pi, sigma, nv, pv, aug, gen);
        drmu(t[j], xi, di, pi, sigma, gen);
        fit(t[j], xi, p, n, x, ftemp);
        for (size_t k = 0; k < n; k++)
            allfit[k] += ftemp[k]; // puts the current tree's prediction back in
    }
    if (dartOn)
    {
        draw_s(nv, lpv, theta, gen);
        draw_theta0(const_theta, theta, lpv, a, b, rho, gen);
        for (size_t j = 0; j < p; j++)
            pv[j] = ::exp(lpv[j]);
    }
}
//--------------------------------------------------
// public functions
void bart::pr() // print to screen
{
    cout << "*****bart object:\n";
    cout << "m: " << m << std::endl;
    cout << "t[0]:\n " << t[0] << std::endl;
    cout << "t[m-1]:\n " << t[m - 1] << std::endl;
    cout << "prior and mcmc info:\n";
    pi.pr();
    if (dart)
    {
        cout << "*****dart prior (On):\n";
        cout << "a: " << a << std::endl;
        cout << "b: " << b << std::endl;
        cout << "rho: " << rho << std::endl;
        cout << "augmentation: " << aug << std::endl;
    }
    else
        cout << "*****dart prior (Off):\n";
    if (p)
        cout << "data set: n,p: " << n << ", " << p << std::endl;
    else
        cout << "data not set\n";
}
