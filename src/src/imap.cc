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

int32_t compute_coverage(const imap_t &imap, int32_t x, int32_t y)
{
	assert(x < y);

	ROI rx(x, x + 1);
	ROI ry(y - 1, y);

	imap_t::const_iterator xi = imap.upper_bound(ROI(x - 1, x));

	if(xi == imap.end()) return 0;

	int32_t s = 0;
	while(upper(xi->first) <= y)
	{
		//printf(" current = [%d, %d), y = %d\n", lower(xi->first), upper(xi->first), y);
		s += upper(xi->first) - lower(xi->first);
		xi++;
		if(xi == imap.end()) break;
	}
	return s;
}

int test_imap()
{
	imap_t imap;
	imap += make_pair(ROI(6, 7), 3);
	imap += make_pair(ROI(1, 3), 3);
	imap += make_pair(ROI(1, 2), 1);
	imap += make_pair(ROI(2, 5), 2);

	imap_t::const_iterator it;
	
	for(it = imap.begin(); it != imap.end(); it++)
	{
		printf("interval: [%d,%d) -> %d\n", lower(it->first), upper(it->first), it->second);
	}

	for(int i = 1; i <= 7; i++)
	{
		it = imap.find(i);
		if(it == imap.end())
		{
			printf("find %d: does not exist\n", i);
		}
		else
		{
			printf("find %d: [%d,%d) -> %d\n", i, lower(it->first), upper(it->first), it->second);
		}
	}

	for(int i = 1; i <= 7; i++)
	{
		it = imap.lower_bound(ROI(i, i + 1));

		if(it == imap.end())
		{
			printf("lower bound %d: does not exist\n", i);
		}
		else
		{
			printf("lower bound %d: [%d,%d) -> %d\n", i, lower(it->first), upper(it->first), it->second);
		}
	}

	for(int i = 1; i <= 7; i++)
	{
		it = imap.upper_bound(ROI(i, i + 1));

		if(it == imap.end())
		{
			printf("upper bound %d: does not exist\n", i);
		}
		else
		{
			printf("upper bound %d: [%d,%d) -> %d\n", i, lower(it->first), upper(it->first), it->second);
		}
	}

	for(int i = 1; i <= 7; i++)
	{
		for(int j = i + 1; j <= 7; j++)
		{
			int s = compute_coverage(imap, i, j);
			printf("coverage [%d,%d) = %d\n", i, j, s);
		}
	}

	return 0;
}
