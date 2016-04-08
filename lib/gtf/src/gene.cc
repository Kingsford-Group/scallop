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

