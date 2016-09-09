#include "slope.h"
#include "util.h"

slope::slope(int t, int x, int s, double sg)
	: type(t), xbin(x), score(s), sigma(sg)
{
}

int slope::print(int index) const
{
	printf("slope %d: type = %s, xbins = %d, score = %d, sigma = %.2lf\n", 
			index, type == 0 ? "5end" : "3end", xbin, score, sigma);
	return 0;
}

bool compare_slope_score(const slope &x, const slope &y)
{
	if(x.score > y.score) return true;
	else return false;
}

