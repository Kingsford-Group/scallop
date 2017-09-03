/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "edge_base.h"
#include <cstdio>

using namespace std;

edge_base::edge_base(int _s, int _t)
	:s(_s), t(_t)
{}

int edge_base::move(int x, int y)
{
	s = x;
	t = y;
	return 0;
}

int edge_base::swap()
{
	int x = s;
	s = t;
	t = x;
	return 0;
}

int edge_base::source() const
{
	return s;
}

int edge_base::target() const
{
	return t;
}

int edge_base::print() const
{
	printf("edge %d -> %d\n", s, t);
	return 0;
}
