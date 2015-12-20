#include <cassert>
#include <cstdio>

#include "bundle.h"

bundle::bundle()
{
	tid = -1;
	lpos = INT32_MAX;
	rpos = 0;
}

int bundle::add_hit(const bam1_core_t &p)
{
	hits.push_back(p);
	int32_t l = p.pos;
	int32_t r = p.pos + p.l_qseq;
	assert(r >= l);

	if(tid == -1) tid = p.tid;
	else assert(tid == p.tid);

	if(l < lpos) lpos = l;
	if(r > rpos) rpos = r;
	return 0;
}

int bundle::clear()
{
	hits.clear();
	tid = -1;
	lpos = INT32_MAX;
	rpos = 0;
	return 0;
}

int bundle::print()
{
	printf("tid = %4d, #hits = %7ld, range = (%8d,%8d)\n", tid, hits.size(), lpos, rpos);
	return 0;
}
