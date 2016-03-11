#include "edge_base.h"
#include <cstdio>

using namespace std;

edge_b::edge_b(int _s, int _t)
	:s(_s), t(_t)
{}

int edge_b::source() const
{
	return s;
}

int edge_b::target() const
{
	return t;
}

int edge_b::print() const
{
	printf("edge %d -> %d\n", s, t);
	return 0;
}
