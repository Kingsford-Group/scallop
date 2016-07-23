#include "vertex_info.h"

vertex_info::vertex_info()
{
	stddev = 1.0;
	length = 0;
	scalor1 = 1;
	scalor2 = 1;
	length1 = 0;
	length2 = 0;
}

vertex_info::vertex_info(int l)
	: length(l)
{
	stddev = 1.0;
	scalor1 = 1;
	scalor2 = 1;
	length1 = 0;
	length2 = 0;
}

vertex_info::vertex_info(const vertex_info &vi)
{
	stddev = vi.stddev;
	length = vi.length;
	scalor1 = vi.scalor1;
	scalor2 = vi.scalor2;
	length1 = vi.length1;
	length2 = vi.length2;
}
