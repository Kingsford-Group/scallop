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
		cout << "add genome " << line.c_str() << endl << flush;
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
		cout << "union genome " << i << endl << flush;
	}
	return 0;
}

int gtfmerge::build_pairwise_intersection(genome1 &gm)
{
	for(int i = 0; i < genomes.size(); i++)
	{
		genomes[i].add_suffix("u" + to_string(i));
	}

	gm.clear();
	for(int i = 0; i < genomes.size(); i++)
	{
		genome1 gx;
		for(int j = i + 1; j < genomes.size(); j++)
		{
			genome1 g1;
			genomes[i].build_intersection(genomes[j], g1);
			gx.build_union(g1);
			printf("union the intersection betwwen %d and %d: (%lu, %lu, %lu)\n", 
					i, j, g1.transcripts.size(), gx.transcripts.size(), gm.transcripts.size());
			fflush(stdout);
		}
		gm.build_union(gx);
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
