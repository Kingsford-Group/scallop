#ifndef __EQUATION_H__
#define __EQUATION_H__

#include <vector>

using namespace std;

class equation
{
public:
	equation();
	equation(double);
	equation(const vector<int> &, const vector<int> &);
	equation(const vector<int> &, const vector<int> &, double);

public:
	int print(int index);

public:
	vector<int> s;		// subs
	vector<int> t;		// subt
	double e;			// erro

	int f;				// 3: resolve vertex 2: fully, 1: partly, 0: none
	int a;				// # adjacent merges
	int d;				// # distant merges
	int w;				// weight
};

bool equation_cmp1(const equation &x, const equation &y);
bool equation_cmp2(const equation &x, const equation &y);
bool equation_cmp3(const equation &x, const equation &y);

#endif
