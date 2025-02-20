

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "heterbart.h"
#include "kevin.h"

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

// TODO: Add multivariate versions of the functions above

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
