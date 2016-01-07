#include <cassert>
#include <cstdio>

#include "bbase.h"

bbase::bbase()
{
	tid = -1;
	chrm = "";
	lpos = INT32_MAX;
	rpos = 0;
}

bbase::~bbase()
{}

int bbase::add_hit(bam_hdr_t *h, bam1_t *b)
{
	// create and store new hit
	hit ht(b);
	hits.push_back(ht);

	// calcuate the boundaries on reference
	if(ht.pos < lpos) lpos = ht.pos;
	if(ht.rpos > rpos) rpos = ht.rpos;

	// set chromsome ID and name
	if(tid == -1)
	{
		tid = ht.tid;
		assert(chrm == "");
		char buf[1024];
		strcpy(buf, h->target_name[ht.tid]);
		chrm = string(buf);
	}
	assert(tid == ht.tid);

	return 0;
}

int bbase::clear()
{
	tid = -1;
	chrm = "";
	lpos = INT32_MAX;
	rpos = 0;
	hits.clear();
	return 0;
}

int32_t bbase::get_tid()
{
	return tid;
}

int32_t bbase::get_rpos()
{
	return rpos;
}

size_t bbase::get_num_hits()
{
	return hits.size();
}
