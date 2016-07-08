#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>

#include "genome.h"
#include "simulator.h"
#include "util.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc != 6)
	{
		cout<<"usage: "<<argv[0]<< " <num-genes> <num_exons> <num_transcripts> <max_expression> <output_file>"<<endl;
		return 0;
	}

	int num_genes = atoi(argv[1]);
	int num_exons = atoi(argv[2]);
	int num_transcripts = atoi(argv[3]);
	int max_expression = atoi(argv[4]);

	srand(time(0));

	genome gm;
	simulator s(num_exons, num_transcripts, max_expression);
	for(int i = 0; i < num_genes; i++)
	{
		string gid = "gene" + tostring(i + 1);
		gene g;
		s.simulate_gene(gid, g);
		gm.add_gene(g);
	}
	gm.write(argv[5]);
    return 0;
}
