#include <cstdio>
#include <cassert>
#include <sstream>
#include <map>

#include "genome.h"

genome::genome(const string &file)
{
	read(file);
	build_index();
}

genome::~genome()
{
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
	map<string, int> m;
	while(fin.getline(line, 102400, '\n'))
	{
		exon ge(line);
		if(ge.feature != "exon") continue;
		if(m.find(ge.gene_id) == m.end())
		{
			gene gg;
			gg.add_exon(ge);
			genes.push_back(gg);
			m.insert(pair<string, int>(ge.gene_id, genes.size() - 1));
		}
		else
		{
			genes[m[ge.gene_id]].add_exon(ge);
		}
	}

	for(int i = 0; i < genes.size(); i++)
	{
		genes[i].build_transcripts();
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

int genome::build_index()
{
	s2i.clear();
	for(int i = 0; i < genes.size(); i++)
	{
		s2i.insert(pair<string, int>(genes[i].gene_id, i));
	}
	return 0;
}

const gene* genome::get_gene(string name) const
{
	map<string, int>::const_iterator it = s2i.find(name);
	if(it == s2i.end()) return NULL;
	int k = it->second;
	return &(genes[k]);
}

