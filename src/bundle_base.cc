/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

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
	strand = '.';
	num_long_reads = 0;
}

bundle_base::~bundle_base()
{}

int bundle_base::add_hit(const hit &ht)
{
	if(ht.is_long_read == true) num_long_reads++;

	// store new hit
	hits.push_back(ht);

	// calcuate the boundaries on reference
	if(ht.pos < lpos) lpos = ht.pos;
	if(ht.rpos > rpos) rpos = ht.rpos;

	// set tid
	if(tid == -1) tid = ht.tid;
	assert(tid == ht.tid);

	// set strand
	if(hits.size() <= 1) strand = ht.strand;
	assert(strand == ht.strand);

	// add matching regions
	int32_t p1 = ht.pos;
	int32_t p2 = -1;
	for(int i = 0; i < ht.spos.size(); i++)
	{
		p2 = high32(ht.spos[i]);
		mmap += make_pair(ROI(p1, p2), 1);
		p1 = low32(ht.spos[i]);
	}
	p2 = ht.rpos;
	mmap += make_pair(ROI(p1, p2), 1);

	/*
	// add intervals
	vector<int64_t> vm;
	vector<int64_t> vi;
	vector<int64_t> vd;
	ht.get_mid_intervals(vm, vi, vd);

	//ht.print();
	for(int k = 0; k < vm.size(); k++)
	{
		int32_t s = high32(vm[k]);
		int32_t t = low32(vm[k]);
		//printf(" add interval %d-%d\n", s, t);
		mmap += make_pair(ROI(s, t), 1);
	}

	for(int k = 0; k < vi.size(); k++)
	{
		int32_t s = high32(vi[k]);
		int32_t t = low32(vi[k]);
		imap += make_pair(ROI(s, t), 1);
	}

	for(int k = 0; k < vd.size(); k++)
	{
		int32_t s = high32(vd[k]);
		int32_t t = low32(vd[k]);
		imap += make_pair(ROI(s, t), 1);
	}
	*/

	return 0;
}

bool bundle_base::overlap(const hit &ht) const
{
	if(mmap.find(ROI(ht.pos, ht.pos + 1)) != mmap.end()) return true;
	if(mmap.find(ROI(ht.rpos - 1, ht.rpos)) != mmap.end()) return true;
	return false;
}

int bundle_base::clear()
{
	tid = -1;
	chrm = "";
	lpos = INT32_MAX;
	rpos = 0;
	strand = '.';
	hits.clear();
	mmap.clear();
	imap.clear();
	num_long_reads = 0;
	return 0;
}
