#include "gtfmerge.h"
#include <string>

int gtfmerge::add_genome(const string &file)
{
	int n = genomes.size();
	genome1 gm;
	genomes.push_back(gm);
	genomes[n].build(file);
	return 0;
}

int gtfmerge::add_genomes(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail()) return 0;

	string line;
	while(getline(fin, line))
	{
		if(line == "") continue;
		add_genome(line);
	}
	return 0;
}

int gtfmerge::build_union(genome1 &gm)
{
	gm.clear();
	for(int i = 0; i < genomes.size(); i++)
	{
		genomes[i].add_suffix("u" + to_string(i));
		gm.build_union(genomes[i]);
	}
	return 0;
}

int gtfmerge::build_pairwise_intersection(genome1 &gm)
{
	for(int i = 0; i < genomes.size(); i++)
	{
		genomes[i].add_suffix("u" + to_string(i));
	}

	vector<genome1> vv;
	for(int i = 0; i < genomes.size(); i++)
	{
		for(int j = i + 1; j < genomes.size(); j++)
		{
			genome1 g1;
			genomes[i].build_intersection(genomes[j], g1);
			vv.push_back(g1);
		}
	}

	gm.clear();
	for(int i = 0; i < vv.size(); i++)
	{
		gm.build_union(vv[i]);
	}
	return 0;
}

int gtfmerge::print()
{
	for(int i = 0; i < genomes.size(); i++)
	{
		genomes[i].print(i);
	}
	return 0;
}
