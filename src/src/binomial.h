#ifndef __BINOMIAL_H__
#define __BINOMIAL_H__

// boost::binomial distribution
#include "boost/math/distributions/binomial.hpp"
using namespace boost::math;

// return the score(transformed from probability)
// that >= x is observed
uint32_t compute_binomial_score(int n, double pr, int x);

#endif
