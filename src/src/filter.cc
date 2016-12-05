#include "filter.h"
#include "config.h"
#include <cassert>
#include <algorithm>

filter::filter(const vector<transcript> &v)
	:trs(v)
{}

int filter::join()
{
	while(true)
	{
		bool b = join_transcripts();
		if(b == false) break;
	}
	return 0;
}

bool filter::join_transcripts()
{
	sort(trs.begin(), trs.end(), transcript_cmp);
	int32_t mind = max_joining_gap + 1;
	int ki = -1, kj = -1;
	for(int i = 0; i < trs.size(); i++)
	{
		int j = locate_next_transcript(i);
		if(j == -1) continue;
		if(trs[i].exons.size() >= 2 && trs[j].exons.size() >= 2) continue;
		int32_t d = trs[j].get_bounds().first - trs[i].get_bounds().second;
		if(d > mind) continue;
		mind = d;
		ki = i;
		kj = j;
	}
	if(ki == -1 || kj == -1) return false;
	if(mind > max_joining_gap) return false;

	if(trs[ki].exons.size() >= 2)
	{
		assert(trs[kj].exons.size() == 1);
		int32_t p1 = trs[ki].get_bounds().second;
		int32_t p2 = trs[kj].get_bounds().second;
		trs[ki].add_exon(p1, p2);
		trs[ki].shrink();
		trs.erase(trs.begin() + kj);
		return true;
	}
	else if(trs[kj].exons.size() >= 2)
	{
		assert(trs[ki].exons.size() == 1);
		int32_t p1 = trs[ki].get_bounds().first;
		int32_t p2 = trs[kj].get_bounds().first;
		trs[kj].add_exon(p1, p2);
		trs[kj].shrink();
		trs.erase(trs.begin() + ki);
		return true;
	}
	else
	{
		assert(trs[ki].exons.size() == 1);
		assert(trs[kj].exons.size() == 1);
		int32_t p1 = trs[ki].get_bounds().first;
		int32_t p2 = trs[kj].get_bounds().first;
		trs[kj].add_exon(p1, p2);
		trs[kj].shrink();
		double cov = 0;
		cov += trs[ki].coverage * trs[ki].length();
		cov += trs[kj].coverage * trs[kj].length();
		cov /= (trs[ki].length() + trs[kj].length());
		trs[kj].coverage = cov;
		trs.erase(trs.begin() + ki);
		return true;
	}

	return true;
}

int filter::locate_next_transcript(int t)
{
	if(t < 0 || t >= trs.size()) return -1;
	PI32 p = trs[t].get_bounds();
	int a = 0;
	int b = trs.size() - 1;
	if(trs[b].get_bounds().first < p.second) return -1;
	while(true)
	{
		assert(a <= b);
		if(a == b) return a;
		int k = (a + b) / 2;
		if(trs[k].get_bounds().first == p.second) return k;
		if(trs[k].get_bounds().first < p.second) a = k + 1;
		if(trs[k].get_bounds().first > p.second) b = k;
	}
	assert(false);
	return -1;
}

bool transcript_cmp(const transcript &x, const transcript &y)
{
	if(x.exons[0].first < y.exons[0].first) return true;
	else return false;
}
