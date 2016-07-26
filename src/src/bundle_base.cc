#include <cassert>
#include <cstdio>
#include <cmath>

#include "bundle_base.h"

bundle_base::bundle_base()
{
	tid = -1;
	chrm = "";
	lpos = INT32_MAX;
	rpos = 0;
	ave_isize = 0;
	strand = '.';
}

bundle_base::~bundle_base()
{}

int bundle_base::add_hit(const hit &ht)
{
	// store new hit
	hits.push_back(ht);

	// calcuate the boundaries on reference
	if(ht.pos < lpos) lpos = ht.pos;
	if(ht.rpos > rpos) rpos = ht.rpos;

	/*
	if(tid == -1)
	{
		tid = ht.tid;
		assert(chrm == "");
		char buf[1024];
		strcpy(buf, h->target_name[ht.tid]);
		chrm = string(buf);
	}
	*/

	// set tid
	if(tid == -1) tid = ht.tid;
	assert(tid == ht.tid);

	// set strand
	if(strand == '.') strand = ht.xs;
	assert(strand == ht.xs);

	if(ht.isize > 0 && ht.isize <= 400)
	{
		double isize = ht.mpos - ht.rpos;
		ave_isize = (ave_isize * (hits.size() - 1) + isize) * 1.0 / hits.size();
	}
	return 0;
}

int bundle_base::clear()
{
	tid = -1;
	chrm = "";
	lpos = INT32_MAX;
	rpos = 0;
	hits.clear();
	ave_isize = 0;
	strand = '.';
	return 0;
}

int bundle_base::set_chrm(const string &s)
{
	chrm = s;
	return 0;
}

int32_t bundle_base::get_tid()
{
	return tid;
}

int32_t bundle_base::get_rpos()
{
	return rpos;
}

size_t bundle_base::get_num_hits()
{
	return hits.size();
}
