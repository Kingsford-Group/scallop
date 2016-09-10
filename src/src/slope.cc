#include "slope.h"
#include "util.h"

slope::slope(int t, int x, int s, double sg, int32_t l, int32_t r)
	: type(t), xbin(x), score(s), sigma(sg), lpos(l), rpos(r)
{
}

int slope::print(int index) const
{
	printf("slope %d: type = %s, xbins = %d, positions = %d-%d, score = %d (%d, %d), sigma = %.2lf (%.2lf, %.2lf)\n", 
			index, type == 0 ? "5end" : "3end", xbin, lpos, rpos, score, score1, score2, sigma, sigma1, sigma2);
	return 0;
}

bool compare_slope_score(const slope &x, const slope &y)
{
	if(x.score > y.score) return true;
	else return false;
}

