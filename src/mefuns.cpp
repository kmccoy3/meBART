

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "heterbart.h"
#include "mefuns.h"
#include <RcppArmadillo.h>


arn gen;

double rnorm(double mu, double sigma)
{
    return mu + sigma * gen.normal();
}

double dnorm(double x, double mu, double sigma)
{
    // TODO: check if this is how the package normally does this
    double tmp = (x - mu) / sigma;
    tmp = -0.5 * tmp * tmp;
    tmp = exp(tmp);
    tmp = tmp / (sigma * RTPI);
    return tmp;
}


double min(double a, double b)
{
    return a < b ? a : b;
}

void MH_ratio(double x)
{
    printf("x: %f\n", x);
    double tmp = exp(x);
    printf("e^x: %f\n", tmp);
    tmp = log(tmp);
    printf("log(e^x): %f\n", tmp);
    return;
}

arma::mat rmvnorm(int n, arma::vec mu, arma::mat sigma) {
    int ncols = sigma.n_cols;
    arma::mat Y = arma::randn(n, ncols);
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

double dmvnorm(arma::vec x, arma::vec mu, arma::mat sigma) {
    
    int dim = x.n_elem;
    double tmp = 1.0;

    tmp /=  pow(RTPI, dim);
    tmp /= sqrt(arma::det(sigma));
    tmp *= exp(-0.5 * arma::as_scalar((x - mu).t() * arma::inv_sympd(sigma) * (x - mu)));

    return tmp;
}