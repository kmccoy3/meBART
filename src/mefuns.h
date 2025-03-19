
#ifndef GUARD_mefuns_h
#define GUARD_mefuns_h

arma::mat rmvnorm(int n, arma::vec mu, arma::mat sigma);
double dmvnorm(arma::vec x, arma::vec mu, arma::mat sigma);

double rnorm(double mu, double sigma);

void MH_ratio(double x);

double dnorm(double x, double mu, double sigma);

double min(double a, double b);

size_t discrete_uniform(const size_t n);

#endif