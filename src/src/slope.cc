#include "slope.h"
#include "util.h"

slope::slope(int t, int l, int r, int s)
	: type(t), lbin(l), rbin(r), score(s)
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
	printf("slope %d: flag = %d, type = %d, bins = (%d - %d), pos = (%d - %d), score = %d, ave = %.2lf, dev = %.2lf\n", 
			index, flag, type, lbin, rbin, lpos, rpos, score, ave, dev);
	return 0;
}

bool compare_slope_score(const slope &x, const slope &y)
{
	if(x.flag < y.flag) return true;
	if(x.flag > y.flag) return false;
	if(x.score > y.score) return true;
	return false;
}

bool compare_slope_pos(const slope &x, const slope &y)
{
	if(x.lpos < y.lpos) return true;
	else return false;
}
