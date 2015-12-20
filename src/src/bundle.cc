#include <cassert>

#include "bundle.h"

bundle::bundle()
{
	lpos = INT32_MAX;
	rpos = 0;
}

int bundle::add_hit(const bam1_core_t &p)
{
	hits.push_back(p);
	int32_t l = p.pos;
	int32_t r = p.pos + p.l_qseq;
	assert(r >= l);

	if(l < lpos) lpos = l;
	if(r > rpos) rpos = r;
	return 0;
}

int bundle::clear()
{
	hits.clear();
	lpos = INT32_MAX;
	rpos = 0;
	return 0;
}

size_t bundle::size()
{
	return hits.size();
}
