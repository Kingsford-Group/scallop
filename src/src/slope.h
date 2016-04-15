#ifndef __SLOPE_H__
#define __SLOPE_H__

class slope
{
public:
	slope(int t, int l, int r, int s);

public:
	int type;			// 5end or 3 end
	int flag;			// margin or middle
	int lbin;
	int rbin;
	int lpos;
	int rpos;
	int score;			// the likelihood of being a slope
	double ave;			// mean of abundance
	double dev;			// dev of abundance

public:
	int distance(const slope &s) const;
	int print(int index) const;
};

bool compare_slope_score(const slope &x, const slope &y);
bool compare_slope_pos(const slope &x, const slope &y);

#endif
