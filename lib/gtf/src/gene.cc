#include <map>
#include <algorithm>
#include <cassert>
#include "gene.h"

int gene::add_transcript(const transcript &t)
{
	transcripts.push_back(t);
	return 0;
}

int gene::clear()
{
	transcripts.clear();
	return 0;
}

int gene::build(const vector<exon> &v)
{
	map<string, int> m;
	for(int i = 0; i < v.size(); i++)
	{
		const exon &e = v[i];
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
	return 0;
}

string gene::get_seqname() const
{
	if(transcripts.size() == 0) return "";
	else return transcripts[0].seqname;
}

string gene::get_gene_id() const
{
	if(transcripts.size() == 0) return "";
	else return transcripts[0].gene_id;
}

int gene::set_gene_id(const string &id) 
{
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcripts[i].gene_id = id;
	}
	return 0;
}

int gene::sort()
{
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcripts[i].sort();
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
