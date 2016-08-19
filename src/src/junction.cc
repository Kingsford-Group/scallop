#include <cstdio>
#include "junction.h"
#include "config.h"

junction::junction()
{}

junction::junction(int64_t _p)
{
	lpos = high32(_p);
	rpos = low32(_p);
	count = 0;
	lexon = -1;
	rexon = -1;
}

junction::junction(int64_t _p, int _c)
{
	lpos = high32(_p);
	rpos = low32(_p);
	count = _c;
	lexon = -1;
	rexon = -1;
}

junction::junction(const junction &sp)
{
	lpos = sp.lpos;
	rpos = sp.rpos;
	count = sp.count;
	lexon = sp.lexon;
	rexon = sp.rexon;
}

bool junction::operator<(const junction &x) const
{
	if(lpos <= x.lpos) return true;
	else return false;
}

int junction::print(int index) const
{
	printf("junction %d: region = [%d, %d), %d -> %d, length = %d, count = %d\n", 
			index, lpos, rpos, lexon, rexon, rpos - lpos, count);
	return 0;
}
