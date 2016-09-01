#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <ctime>
#include <cassert>
#include <sstream>

#include "assembler.h"
#include "config.h"
#include "bundle.h"
#include "scallop1.h"
#include "gtf.h"
#include "genome.h"
#include "subsetsum4.h"

using namespace std;

int simulate();
int process_bam();
int process_gtf();
int process_sgr();

int main(int argc, const char **argv)
{
	srand(time(0));
	parse_arguments(argc, argv);
	print_parameters();

	if(input_file == "") simulate();

	string s = input_file.substr(input_file.size() - 3, 3);
	if(s == "bam") process_bam();
	if(s == "gtf") process_gtf();
	if(s == "sgr") process_sgr();

	return 0;

}

int simulate()
{
	if(output_file == "") return 0;
	splice_graph sg;
	sg.simulate(simulation_num_vertices, simulation_num_edges, simulation_max_edge_weight);
	sg.write(output_file);
	return 0;
}

int process_bam()
{
	assembler asmb;
	asmb.assemble();
	return 0;
}

int process_gtf()
{
	genome g(input_file);

	ofstream fout;
	if(output_file != "") fout.open(output_file);

	for(int i = 0; i < g.genes.size(); i++)
	{
		gtf gg(g.genes[i]);

		if(gg.transcripts.size() < min_gtf_transcripts_num) continue;

		string name = gg.get_gene_id();

		if(fixed_gene_name != "" && name != fixed_gene_name) continue;

		splice_graph gr;
		gg.build_splice_graph(gr);

		scallop1 sc(name, gr);
		sc.assemble();

		if(output_file == "") continue;
		gg.output_gtf(fout, sc.paths, algo);
	}

	if(output_file != "") fout.close();

	return 0;
}

int process_sgr()
{
	splice_graph sg;
	sg.build(input_file);

	scallop1 sc("shao", sg);
	sc.assemble();

	return 0;
}
