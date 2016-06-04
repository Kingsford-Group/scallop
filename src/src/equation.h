#ifndef __EQUATION_H__
#define __EQUATION_H__

#include <vector>

using namespace std;

class equation
{
public:
	equation(const vector<int> &, const vector<int> &);
	equation(const vector<int> &, const vector<int> &, double);

	bool operator< (const equation &e) const;

public:
	vector<int> s;
	vector<int> t;
	double e;
};

#endif
