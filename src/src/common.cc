#include "common.h"

#include <boost/math/distributions/binomial.hpp>

using namespace boost::math;

int compute_overlap(const imap_t &imap, int32_t p)
{
	imap_t::const_iterator it = imap.find(p);
	if(it == imap.end()) return 0;
	return it->second;
}

int cumulate_overlap(const imap_t &imap, int32_t x, int32_t y)
{
	int s = 0;
	for(int k = x; k < y; k++)
	{
		s += compute_overlap(imap, k);
	}
	return s;
}

uint32_t compute_binomial_score(int n, double pr, int x)
{
	// compute the score that observed >= x among n trials
	assert(x >= 0 && x <= n);
	if(x == 0) return UINT32_MAX;

	binomial_distribution<> b(n, pr);
	double p = cdf(complement(b, x - 1));
	return (uint32_t)(-10.0 * log10(p));
}
