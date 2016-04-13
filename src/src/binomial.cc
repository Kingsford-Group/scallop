#include "binomial.h"

uint32_t compute_binomial_score(int n, double pr, int x)
{
	// compute the score that observed >= x among n trials
	assert(x >= 0 && x <= n);
	if(x == 0) return UINT32_MAX;

	binomial_distribution<> b(n, pr);
	double p = cdf(complement(b, x - 1));
	return (uint32_t)(-100.0 * log10(p));
}
