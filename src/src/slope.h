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
	bool operator < (const slope &s) const;

public:
	int distance(const slope &s) const;

};

#endif
