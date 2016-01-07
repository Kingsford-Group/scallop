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
