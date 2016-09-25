#include "config.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

//// parameters
// for bam file and reads
int min_flank_length = 5;
int32_t min_bundle_gap = 50;
int min_num_hits_in_bundle = 20;
uint32_t min_mapping_quality = 1;
int32_t min_splice_boundary_hits = 1;

// for identifying subgraphs
double max_indel_ratio = 0.2;
int32_t min_subregion_gap = 3;
int32_t min_subregion_length = 15;
double min_subregion_overlap = 2;

// for splice graph
double min_edge_weight = 2.9;
double max_ignorable_edge_weight = 20.0;
double min_transcript_coverage = 10.0;
double min_splice_graph_coverage = 20.0;
double max_split_error_ratio = 0.15;

// for identifying new boundaries
int min_boundary_length = 80;
int min_boundary_score = 1000;
double min_boundary_sigma = 4.0;

// for subsetsum and router
int max_dp_table_size = 10000;
int min_router_count = 1;

// for simulation
int simulation_num_vertices = 0;
int simulation_num_edges = 0;
int simulation_max_edge_weight = 0;

// input and output
string algo = "shao";
string input_file;
string ref_file;
string ref_file1;
string ref_file2;
string output_file;

// for controling
int max_num_bundles = -1;
int32_t average_read_length = 100;
bool strand_reverse = false;
bool ignore_single_exon_transcripts = false;
bool output_tex_files = false;
string fixed_gene_name = "";
int min_gtf_transcripts_num = 0;
bool fast_mode = true;

int print_parameters()
{
	printf("parameters:\n");

	// for bam file and reads
	printf("min_flank_length = %d\n", min_flank_length);
	printf("min_bundle_gap = %d\n", min_bundle_gap);
	printf("min_num_hits_in_bundle = %d\n", min_num_hits_in_bundle);
	printf("min_mapping_quality = %d\n", min_mapping_quality);
	printf("min_splice_boundary_hits = %d\n", min_splice_boundary_hits);

	// for identifying subgraphs
	printf("max_indel_ratio = %.2lf\n", max_indel_ratio);
	printf("min_subregion_gap = %d\n", min_subregion_gap);
	printf("min_subregion_length = %d\n", min_subregion_length);
	printf("min_subregion_overlap = %.2lf\n", min_subregion_overlap);

	// for splice graph
	printf("min_edge_weight = %.2lf\n", min_edge_weight);
	printf("max_ignorable_edge_weight = %.2lf\n", max_ignorable_edge_weight);
	printf("min_transcript_coverage = %.2lf\n", min_transcript_coverage);
	printf("min_splice_graph_coverage = %.2lf\n", min_splice_graph_coverage);
	printf("max_split_error_ratio = %.2lf\n", max_split_error_ratio);

	// for identifying new boundaries
	printf("min_boundary_length = %d\n", min_boundary_length);
	printf("min_boundary_score = %d\n", min_boundary_score);
	printf("min_boundary_signma = %.2lf\n", min_boundary_sigma);

	// for subsetsum and router
	printf("max_dp_table_size = %d\n", max_dp_table_size);
	printf("min_router_count = %d\n", min_router_count);

	// for simulation
	printf("simulation_num_vertices = %d\n", simulation_num_vertices);
	printf("simulation_num_edges = %d\n", simulation_num_edges);
	printf("simulation_max_edge_weight = %d\n", simulation_max_edge_weight);

	// for input and output
	printf("algo = %s\n", algo.c_str());
	printf("input_file = %s\n", input_file.c_str());
	printf("ref_file = %s\n", ref_file.c_str());
	printf("ref_file1 = %s\n", ref_file1.c_str());
	printf("ref_file2 = %s\n", ref_file2.c_str());
	printf("output_file = %s\n", output_file.c_str());

	// for controling
	printf("max_num_bundles = %d\n", max_num_bundles);
	printf("average_read_length = %d\n", average_read_length);
	printf("strand_reverse = %c\n", strand_reverse ? 'T' : 'F');
	printf("ignore_single_exon_transcripts = %c\n", ignore_single_exon_transcripts ? 'T' : 'F');
	printf("output_tex_files = %c\n", output_tex_files ? 'T' : 'F');
	printf("fixed_gene_name = %s\n", fixed_gene_name.c_str());
	printf("min_gtf_transcripts_num = %d\n", min_gtf_transcripts_num);
	printf("fast_mode = %c\n", fast_mode ? 'T' : 'F');

	printf("\n");

	return 0;
}

bool parse_arguments(int argc, const char ** argv)
{
	output_tex_files = false;
	bool b = false;
	for(int i = 1; i < argc; i++)
	{
		if(string(argv[i]) == "-a")
		{
			algo = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-i")
		{
			input_file = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-o")
		{
			output_file = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-r")
		{
			ref_file = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-r1")
		{
			ref_file1 = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-r2")
		{
			ref_file2 = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-g")
		{
			fixed_gene_name = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-t")
		{
			output_tex_files = true;
		}
		else if(string(argv[i]) == "-RF")
		{
			strand_reverse = true;
		}
		else if(string(argv[i]) == "-m")
		{
			ignore_single_exon_transcripts = true;
		}
		else if(string(argv[i]) == "-x")
		{
			min_edge_weight = atof(argv[i + 1]);
			i++;
		}
	}

	return b;
}
