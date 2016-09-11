#include "slope.h"
#include "util.h"

slope::slope(int t, int lb, int rb, int s, double sg)
	: type(t), lbin(lb), rbin(rb), score(s), sigma(sg)
{
	score0 = score1 = score2 = 0;
	sigma0 = sigma1 = sigma2 = 0;
}

int slope::print(int index) const
{
	printf("slope %d: type = %d, bins = %d-%d, positions = %d-%d, score = %d (%d, %d, %d), sigma = %.2lf (%.2lf, %.2lf, %.2lf)\n", 
			index, type, lbin, rbin, lpos, rpos, score, score0, score1, score2, sigma, sigma0, sigma1, sigma2);
	return 0;
}

bool compare_slope_score(const slope &x, const slope &y)
{
	if(x.score > y.score) return true;
	else return false;
}

bool compare_slope_pos(const slope &x, const slope &y)
{
	if(x.lpos < y.lpos) return true;
	else return false;
}
