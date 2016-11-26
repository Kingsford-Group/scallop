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
	assert(g2i.find(g.get_gene_id()) == g2i.end());
	g2i.insert(pair<string, int>(g.get_gene_id(), genes.size()));
	genes.push_back(g);
	return 0;
}

int genome::read(const string &file)
{
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
		if(ge.feature == "transcript")
		{
			if(g2i.find(ge.gene_id) == g2i.end())
			{
				gene gg;
				gg.add_transcript(ge);
				g2i.insert(pair<string, int>(ge.gene_id, genes.size()));
				genes.push_back(gg);
			}
			else
			{
				int k = g2i[ge.gene_id];
				genes[k].add_exon(ge);
			}
		}
		else if(ge.feature != "exon")
		{
			assert(g2i.find(ge.gene_id) != g2i.end());
			int k = g2i[ge.gene_id];
			genes[k].add_exon(ge);
		}
	}

	for(int i = 0; i < genes.size(); i++)
	{
		genes[i].sort();
		genes[i].shrink();
		genes[i].assign_coverage_ratio();
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
