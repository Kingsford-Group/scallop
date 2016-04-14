#ifndef __SLOPE_H__
#define __SLOPE_H__

class slope
{
public:
	slope(int t, int l, int r, int s);

public:
	int type;
	int lpos;
	int rpos;
	int score;
	double dev;

public:
	int distance(const slope &s) const;
	int print(int index) const;
};

bool compare_slope_score(const slope &x, const slope &y);
bool compare_slope_pos(const slope &x, const slope &y);

#endif
