/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>
#include <map>

#include "genome.h"
#include "util.h"

genome::genome()
{}

genome::genome(const string &file)
{
	read(file);
}

genome::~genome()
{
}

int genome::add_gene(const gene &g)
{
	//printf("add gene %s with %lu transcripts\n", g.get_gene_id().c_str(), g.transcripts.size());
	assert(g2i.find(g.get_gene_id()) == g2i.end());
	g2i.insert(pair<string, int>(g.get_gene_id(), genes.size()));
	genes.push_back(g);
	return 0;
}

int genome::read(const string &file)
{
	if(file == "") return 0;

	ifstream fin(file.c_str());
	if(fin.fail())
	{
		printf("open file %s error\n", file.c_str());
		return 0;
	}

	char line[102400];
	
	genes.clear();
	g2i.clear();
	while(fin.getline(line, 102400, '\n'))
	{
		item ge(line);
		if(g2i.find(ge.gene_id) == g2i.end())
		{
			gene gg;
			if(ge.feature == "transcript") gg.add_transcript(ge);
			else if(ge.feature == "exon") gg.add_exon(ge);
			g2i.insert(pair<string, int>(ge.gene_id, genes.size()));
			genes.push_back(gg);
		}
		else
		{
			int k = g2i[ge.gene_id];
			if(ge.feature == "transcript") genes[k].add_transcript(ge);
			else if(ge.feature == "exon") genes[k].add_exon(ge);
		}
	}

	for(int i = 0; i < genes.size(); i++)
	{
		genes[i].sort();
		genes[i].shrink();
	}

	return 0;
}

int genome::write(const string &file) const
{
	ofstream fout(file.c_str());
	for(int i = 0; i < genes.size(); i++)
	{
		genes[i].write(fout);
	}
	fout.close();
	return 0;
}

int genome::sort()
{
	for(int i = 0; i < genes.size(); i++)
	{
		genes[i].sort();
	}
	return 0;
}

const gene* genome::get_gene(string name) const
{
	map<string, int>::const_iterator it = g2i.find(name);
	if(it == g2i.end()) return NULL;
	int k = it->second;
	return &(genes[k]);
}

const gene* genome::locate_gene(const string &chrm, const PI32 &p) const
{
	assert(p.first <= p.second);
	const gene * x = NULL;
	int32_t oo = 0;
	for(int i = 0; i < genes.size(); i++)
	{
		const gene &g = genes[i];
		if(g.get_seqname() != chrm) continue;
		PI32 b = g.get_bounds();
		assert(b.first <= b.second);
		int32_t o = compute_overlap(p, b);
		if(o > 0 && o > oo)
		{
			x = &(genes[i]);
			oo = o;
		}
	}
	return x;
}

int genome::assign_RPKM(double factor)
{
	for(int i = 0; i < genes.size(); i++)
	{
		genes[i].assign_RPKM(factor);
	}
	return 0;
}

int genome::assign_TPM_by_RPKM()
{
	double sum = 0;
	for(int i = 0; i < genes.size(); i++)
	{
		vector<transcript> &v = genes[i].transcripts;
		for(int k = 0; k < v.size(); k++)
		{
			transcript &t = v[k];
			sum += t.RPKM;
		}
	}

	for(int i = 0; i < genes.size(); i++)
	{
		vector<transcript> &v = genes[i].transcripts;
		for(int k = 0; k < v.size(); k++)
		{
			transcript &t = v[k];
			t.TPM = t.RPKM * 1e6 / sum;
		}
	}
	return 0;
}

int genome::assign_TPM_by_FPKM()
{
	double sum = 0;
	for(int i = 0; i < genes.size(); i++)
	{
		vector<transcript> &v = genes[i].transcripts;
		for(int k = 0; k < v.size(); k++)
		{
			transcript &t = v[k];
			sum += t.FPKM;
		}
	}

	for(int i = 0; i < genes.size(); i++)
	{
		vector<transcript> &v = genes[i].transcripts;
		for(int k = 0; k < v.size(); k++)
		{
			transcript &t = v[k];
			t.TPM = t.FPKM * 1e6 / sum;
		}
	}
	return 0;
}

int genome::filter_single_exon_transcripts()
{
	for(int i = 0; i < genes.size(); i++)
	{
		genes[i].filter_single_exon_transcripts();
	}
	return 0;
}

int genome::filter_low_coverage_transcripts(double min_coverage)
{
	for(int i = 0; i < genes.size(); i++)
	{
		genes[i].filter_low_coverage_transcripts(min_coverage);
	}
	return 0;
}

vector<transcript> genome::collect_transcripts() const
{
	vector<transcript> vv;
	for(int i = 0; i < genes.size(); i++)
	{
		const vector<transcript> &v = genes[i].transcripts;
		vv.insert(vv.end(), v.begin(), v.end());
	}
	return vv;
}
