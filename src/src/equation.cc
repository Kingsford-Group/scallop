#include "equation.h"

equation::equation(const vector<int> &_s, const vector<int> &_t)
	: s(_s), t(_t)
{
}

equation::equation(const vector<int> &_s, const vector<int> &_t, double _e)
	: s(_s), t(_t), e(_e)
{}

bool equation::operator< (const equation &eqn) const
{
	if(e < eqn.e - 0.00001) return true;
	else if(e > eqn.e + 0.00001) return false;
	else if(s.size() < eqn.s.size()) return true;
	else return false;
}
