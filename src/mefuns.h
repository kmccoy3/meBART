
#ifndef GUARD_mefuns_h
#define GUARD_mefuns_h

arma::mat rmvnorm(int n, arma::vec mu, arma::mat sigma);

double dmvnorm(arma::vec x, arma::vec mu, arma::mat sigma);

double rnorm(double mu, double sigma);

double dnorm(double x, double mu, double sigma);

RcppExport SEXP cdnorm(SEXP x, SEXP mu, SEXP sigma);

#endif