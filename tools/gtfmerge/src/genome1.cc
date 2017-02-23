#include "genome1.h"
#include <cassert>

genome1::genome1()
{}

genome1::genome1(const string &file)
	: genome(file)
{}

int genome1::build()
{
	collect_multiexon_transcripts();
	build_intron_index();
	return 0;
}

int genome1::collect_multiexon_transcripts()
{
	for(int i = 0; i < genes.size(); i++)
	{
		gene &g = genes[i];
		for(int k = 0; k < g.transcripts.size(); k++)
		{
			transcript &t = g.transcripts[k];
			if(t.exons.size() <= 1) continue;
			transcripts.push_back(t);
		}
	}
	return 0;
}

int genome1::build_intron_index()
{
	intron_index.clear();
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcript &t = transcripts[i];
		assert(t.exons.size() >= 2);
		PI32 p = t.get_first_intron();
		if(intron_index.find(p) == intron_index.end()) 
		{
			set<int> s;
			s.insert(i);
			intron_index.insert(PPIS(p, s));
		}
		else
		{
			intron_index[p].insert(i);
		}
	}
	return 0;
}

int genome1::query(const transcript &t, const set<int> &fb)
{
	if(t.exons.size() <= 1) return -1;
	PI32 p = t.get_first_intron();
	if(intron_index.find(p) == intron_index.end()) return -1;
	set<int> s = intron_index[p];
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		transcript &x = transcripts[k];
		if(x.strand != t.strand) continue;
		if(x.exons.size() != t.exons.size()) continue;
		if(x.intron_chain_match(t) == false) continue;
		if(fb.find(k) != fb.end()) continue;
		return k;
	}
	return -1;
}

int genome1::print(int index)
{
	printf("genome %d: %lu transcripts, %lu distinct first intron\n", index, transcripts.size(), intron_index.size());
	return 0;
}
