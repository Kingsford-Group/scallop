#include "gtfmerge.h"

int gtfmerge::add_genome(const string &file)
{
	int n = genomes.size();
	genome1 gm;
	genomes.push_back(gm);
	genomes[n].build(file);
	return 0;
}

int gtfmerge::build_union(genome1 &gm)
{
	gm.clear();
	for(int i = 0; i < genomes.size(); i++)
	{
		gm.build_union(genomes[i]);
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
