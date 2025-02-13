

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "heterbart.h"
#include "kevin.h"

arn gen;

double kevin_func(double mu, double sigma)
{
    return mu + sigma * gen.normal();
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