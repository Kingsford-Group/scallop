#include "equation.h"
#include <cstdio>

equation::equation(double _e)
	:e(_e)
{}

equation::equation(const vector<int> &_s, const vector<int> &_t)
	: s(_s), t(_t)
{
	e = 0;
	b = false;
	d = 0;
}

equation::equation(const vector<int> &_s, const vector<int> &_t, double _e)
	: s(_s), t(_t), e(_e)
{
	b = false;
	d = 0;
}

int equation::print(int index)
{
	printf("equation %d: (%lu, %lu) edges, error = %.1lf, b = %c, distant = %d. ", 
			index, s.size(), t.size(), e, b ? 'T' : 'F', d);

	printf("S = (%d", s[0]);
	for(int i = 1; i < s.size(); i++) printf(", %d", s[i]);
	printf("), T = (%d", t[0]);
	for(int i = 1; i < t.size(); i++) printf(", %d", t[i]);
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
	if(x.b == true && y.b == false) return true;
	else if(x.b == false && y.b == true) return false;
	else if(x.d < y.d) return true;
	else if(x.d > y.d) return false;
	else if(x.s.size() + x.t.size() < y.s.size() + y.t.size()) return true;
	else return false;
}
