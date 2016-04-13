#include "slope.h"
#include "util.h"

slope::slope(int t, int l, int r, int s)
	: type(t), lpos(l), rpos(r), score(s)
{
	dev = 0;
}

bool slope::operator< (const slope &s) const
{
	if(score > s.score) return true;
	else return false;
}

int slope::distance(const slope &s) const
{
	int o = compute_overlap(PI(lpos, rpos), PI(s.lpos, s.rpos));
	return 0 - o;
}
