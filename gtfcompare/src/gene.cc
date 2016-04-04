#include <map>
#include <algorithm>
#include <cassert>
#include "gene.h"

int gene::add_exon(const exon &ge)
{
	seqname = ge.seqname;
	gene_id = ge.gene_id;
	exons.push_back(ge);
	return 0;
}

int gene::build_transcripts()
{
	map<string, int> m;
	for(int i = 0; i < exons.size(); i++)
	{
		exon &e = exons[i];
		if(m.find(e.transcript_id) == m.end())
		{
			transcript t;
			t.add_exon(e);
			m.insert(pair<string, int>(e.transcript_id, transcripts.size()));
			transcripts.push_back(t);
		}
		else
		{
			transcripts[m[e.transcript_id]].add_exon(e);
		}
	}

	for(int i = 0; i < transcripts.size(); i++)
	{
		transcripts[i].sort();
	}
	return 0;
}

int gene::print() const
{
	for(int i = 0; i < exons.size(); i++)
	{
		exons[i].print();
	}
	return 0;
}

int gene::write(ofstream &fout) const
{
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcripts[i].write(fout);
	}
	return 0;
}

int gene::compare(const gene &g) const
{
	int cnt = 0;
	vector<bool> v;
	v.assign(g.transcripts.size(), false);
	for(int i = 0; i < transcripts.size(); i++)
	{
		for(int j = 0; j < g.transcripts.size(); j++)
		{
			if(v[j] == true) continue;
			if(transcripts[i].equal(g.transcripts[j]) == false) continue;
			v[j] = true;
			cnt++;
			break;
		}
	}
	return cnt;
}
