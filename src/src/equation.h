#ifndef __EQUATION_H__
#define __EQUATION_H__

#include <vector>

using namespace std;

class equation
{
public:
	equation(const vector<int> &, const vector<int> &);
	equation(const vector<int> &, const vector<int> &, double);

public:
	int print(int index);

public:
	vector<int> s;		// subs
	vector<int> t;		// subt
	double e;			// erro

	bool b;				// whether fail or not
	int d;				// # distant merges
};

bool equation_cmp1(const equation &x, const equation &y);
bool equation_cmp2(const equation &x, const equation &y);

#endif
