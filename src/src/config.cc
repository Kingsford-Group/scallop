#include "config.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

//// parameters

// for bam file
int32_t min_bundle_gap = 100;
int min_num_hits_in_bundle = 20;
int32_t min_splice_boundary_hits = 1;
uint32_t min_max_splice_boundary_qual = 3;
double min_average_overlap = 2;
int min_max_region_overlap = 5;
double min_region_coverage = 0.5;
int max_num_bundles = -1;
int slope_bin_size = 10;
int slope_min_bin_num = 9;
int slope_max_bin_num = 32;
int slope_min_distance = 100;
int slope_min_score = 100;
int slope_flexible_size = 20;
double slope_acceptance_dev_decrease = 0.05;
int32_t average_read_length = 100;
int32_t average_slope_length = 240;
int pseudo_length_count = 10;
double min_boundary_edge_weight_ratio = 0.05;

// for algorithm
int max_dp_table_size = 10000;
int max_num_subsetsum_solutions = 10;
double max_equation_error_ratio = 0.2;

// for simulation
int simulation_num_vertices = 0;
int simulation_num_edges = 0;
int simulation_max_edge_weight = 0;

//// from command line
string algo;
string input_file;
string output_file;

bool output_tex_files = false;
string fixed_gene_name = "";
int min_gtf_transcripts_num = 0;
bool fast_mode = true;

int print_parameters()
{
	// for bam file
	printf("parameters:\n");
	printf("min_bundle_gap = %d\n", min_bundle_gap);
	printf("min_num_hits_in_bundle = %d\n", min_num_hits_in_bundle);
	printf("min_splice_boundary_hits = %d\n", min_splice_boundary_hits);
	printf("min_max_splice_boundary_qual = %d\n", min_max_splice_boundary_qual);
	printf("min_average_overlap = %.2lf\n", min_average_overlap);
	printf("min_max_region_overlap = %d\n", min_max_region_overlap);
	printf("min_region_coverage = %.2lf\n", min_region_coverage);
	printf("max_num_bundles = %d\n", max_num_bundles);
	printf("slope_bin_size = %d\n", slope_bin_size);
	printf("slope_min_bin_num = %d\n", slope_min_bin_num);
	printf("slope_min_distance = %d\n", slope_min_distance);
	printf("slope_min_score = %d\n", slope_min_score);
	printf("slope_flexible_size = %d\n", slope_flexible_size);
	printf("slope_acceptance_dev_decrease = %.2lf\n", slope_acceptance_dev_decrease);
	printf("average_read_length = %d\n", average_read_length);
	printf("average_slope_length = %d\n", average_slope_length);
	printf("pseudo_length_count = %d\n", pseudo_length_count);
	printf("min_boundary_edge_weight_ratio = %.2lf\n", min_boundary_edge_weight_ratio);

	// for algorithm
	printf("max_dp_table_size = %d\n", max_dp_table_size);
	printf("max_num_subsetsum_solutions = %d\n", max_num_subsetsum_solutions);
	printf("max_equation_error_ratio = %.2lf\n", max_equation_error_ratio);

	// for simulation
	printf("simulation_num_vertices = %d\n", simulation_num_vertices);
	printf("simulation_num_edges = %d\n", simulation_num_edges);
	printf("simulation_max_edge_weight = %d\n", simulation_max_edge_weight);

	// for command
	printf("algo = %s\n", algo.c_str());
	printf("input_file = %s\n", input_file.c_str());
	printf("output_file = %s\n", output_file.c_str());
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
		else if(string(argv[i]) == "-o")
		{
			output_file = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-i")
		{
			input_file = string(argv[i + 1]);
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
		else if(string(argv[i]) == "-f")
		{
			fast_mode = true;
		}
		else if(string(argv[i]) == "-x")
		{
			slope_min_score = atoi(argv[i + 1]);
			i++;
		}
	}

	return b;
}

