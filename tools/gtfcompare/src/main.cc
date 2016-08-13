#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include "genome.h"
#include "compare.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc != 4)
	{
		cout<<"usage: "<<argv[0]<< " <gtf-file-1> <gtf-file-2> <mode: (1: strict, 2: intron_chain)>"<<endl;
		return 0;
	}

	int mode = atoi(argv[3]);

	genome g1(argv[1]);
	genome g2(argv[2]);

	remove_single_exon_transcripts(g1);
	remove_single_exon_transcripts(g2);

	if(mode == 1) compare_genome1(g1, g2);
	if(mode == 2) compare_genome2(g1, g2);

    return 0;
}
