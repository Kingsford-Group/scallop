#include "vertex_info.h"

vertex_info::vertex_info()
{
	stddev = 1.0;
	length = 0;
	sdist = -1;
	tdist = -1;
	reliability = 0;
	infer = false;
	type = -1;
	adjust = false;
	lpos = 0;
	rpos = 0;
	pos = 0;
}

vertex_info::vertex_info(int l)
	: length(l)
{
	stddev = 1.0;
	sdist = -1;
	tdist = -1;
	reliability = 0;
	infer = false;
	type = -1;
	adjust = false;
	lpos = 0;
	rpos = 0;
	pos = 0;
}

vertex_info::vertex_info(const vertex_info &vi)
{
	stddev = vi.stddev;
	length = vi.length;
	sdist = vi.sdist;
	tdist = vi.tdist;
	reliability = vi.reliability;
	infer = vi.infer;
	type = vi.type;
	adjust = vi.adjust;
	lpos = vi.lpos;
	rpos = vi.rpos;
	pos = vi.pos;
}
