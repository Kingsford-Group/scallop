#include "segment.h"
#include "config.h"
#include "util.h"
#include "binomial.h"
#include <algorithm>

using namespace std;

segment::segment(const split_interval_map *_imap)
	:imap(_imap)
{
}

segment::~segment()
{}

int segment::clear()
{
	pexons.clear();
	bins.clear();
	return 0;
}

int segment::build()
{
	return 0;
}

int segment::add_partial_exon(int k, const partial_exon &pe, double w)
{
	plist.push_back(k);
	pexons.push_back(pe);
	offsets.push_back(w);
	return 0;
}

int segment::print(int index) const
{
	printf("segment %d: \n", index);
	for(int i = 0; i < pexons.size(); i++)
	{
		printf(" index = %d, offset = %.2lf, ", plist[i] + 1, offsets[i]);
		pexons[i].print(i);
	}
	return 0;
}
