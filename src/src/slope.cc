#include "slope.h"
#include "util.h"

slope::slope(int t, int l, int r, int s)
	: type(t), lpos(l), rpos(r), score(s)
{
	dev = 0;
}

int slope::distance(const slope &s) const
{
	int o = compute_overlap(PI(lpos, rpos), PI(s.lpos, s.rpos));
	return 0 - o;
}

int slope::print(int index) const
{
	printf("slope %d: type = %d, (%d - %d), score = %d, dev = %.2lf\n", index, type, lpos, rpos, score, dev);
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

