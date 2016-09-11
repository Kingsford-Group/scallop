#ifndef __SLOPE_H__
#define __SLOPE_H__

#include <stdint.h>

class slope
{
public:
	slope(int t, int lbin, int rbin, int s, double sg);

public:
	int type;			// 5end or 3end
	int lbin;			// index in bsets
	int rbin;			// index in bsets
	int score;			// the likelihood of being a slope
	double sigma;		// variance sigma
	int32_t lpos;		// left position
	int32_t rpos;		// right position
	int score0, score1, score2;
	double sigma0, sigma1, sigma2;

public:
	int print(int index) const;
};

bool compare_slope_score(const slope &x, const slope &y);
bool compare_slope_pos(const slope &x, const slope &y);

#endif
