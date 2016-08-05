#include "equation.h"
#include "util.h"
#include <cstdio>

equation::equation()
{
	e = 0;
	f = 0;
	a = 0;
	d = 0;
	w = 0;
}

equation::equation(double _e)
	:e(_e)
{
	f = 0;
	a = 0;
	d = 0;
	w = 0;
}

equation::equation(const vector<int> &_s, const vector<int> &_t)
	: s(_s), t(_t)
{
	e = 0;
	f = 0;
	d = 0;
	a = 0;
	w = 0;
}

equation::equation(const vector<int> &_s, const vector<int> &_t, double _e)
	: s(_s), t(_t), e(_e)
{
	f = 0;
	d = 0;
	a = 0;
	w = 0;
}

int equation::print(int index)
{
	printf("equation %3d: (%2lu, %2lu) edges, error = %3.1lf, f = %d, w = %d, adjacent = %2d, distant = %2d. ", 
			index, s.size(), t.size(), e, f, w, a, d);

	printf("S = (");
	printv(s);
	printf("), T = (");
	printv(t);
	printf(")\n");
	
	return 0;
}

bool equation_cmp1(const equation &x, const equation &y) 
{
	if(x.e < y.e - 0.00001) return true;
	else if(x.e > y.e + 0.00001) return false;
	else if(x.s.size() + x.t.size() < y.s.size() + y.t.size()) return true;
	else return false;
}

bool equation_cmp2(const equation &x, const equation &y) 
{
	if(x.f > y.f) return true;
	if(x.f < y.f) return false;
	if(x.d < y.d) return true;
	if(x.d > y.d) return false;
	if(x.s.size() + x.t.size() < y.s.size() + y.t.size()) return true;
	if(x.s.size() + x.t.size() > y.s.size() + y.t.size()) return false;
	if(x.w < y.w) return true;
	return false;
}
