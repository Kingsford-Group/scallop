#ifndef __SLOPE_H__
#define __SLOPE_H__

class slope
{
public:
	slope(int t, int x, int s, double sg);

public:
	int type;			// 5end or 3end
	int xbin;			// index in bsets
	int score;			// the likelihood of being a slope
	double sigma;		// variance sigma

public:
	int print(int index) const;
};

bool compare_slope_score(const slope &x, const slope &y);

#endif
