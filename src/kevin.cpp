

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


