#include "imap.h"

int compute_overlap(const imap_t &imap, int32_t p)
{
	imap_t::const_iterator it = imap.find(p);
	if(it == imap.end()) return 0;
	return it->second;
}

int cumulate_overlap(const imap_t &imap, int32_t x, int32_t y, int32_t t)
{
	int s = 0;
	for(int k = x; k < y; k+=t)
	{
		s += compute_overlap(imap, k);
	}
	return s;
}

int maximum_overlap(const imap_t &imap, int32_t x, int32_t y, int32_t t)
{
	int s = 0;
	for(int k = x; k < y; k+=t)
	{
		int z = compute_overlap(imap, k);
		if(z > s) s = z;
	}
	return s;
}

int compute_coverage(const imap_t &imap, int32_t x, int32_t y, int32_t t)
{
	int s = 0;
	for(int k = x; k < y; k+=t)
	{
		if(compute_overlap(imap, k) >= 1) s++;
	}
	return s;
}
