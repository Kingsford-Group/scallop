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
int min_flank_length = 1;
int32_t average_read_length = 100;
int32_t average_slope_length = 240;
int32_t min_bundle_gap = 100;
int32_t min_subregion_gap = 10;
int min_num_hits_in_bundle = 10;
int32_t min_splice_boundary_hits = 1;
uint32_t min_mapping_quality = 5;
double max_indel_ratio = 0.2;
double min_average_overlap = 3.0;
int32_t min_subregion_length = 50;
int min_max_region_overlap = 5;
double min_region_coverage = 0.5;
int max_num_bundles = -1;
int tail_coverage = 8;
bool identify_slopes = false;
int slope_bin_size = 5;
int slope_min_score = 30;
int slope_extend_score = 24;
int slope_min_bin_num = 16;
int slope_flexible_bin_num = 2;
double slope_acceptance_sigma = 2.0;
int pseudo_length_count = 10;
bool use_paired_end = false;
int max_equations_each_iteration = 50;

// for algorithm
double join_min_reliability = 0.6;
double infer_min_reliability = 0.6;
double infer_root_reliability = 0.4;
int max_dp_table_size = 10000;
int max_num_subsetsum_solutions = 10;
double max_equation_error_ratio = 0.1;
double max_gateway_error_ratio = 0.3;
double min_boundary_edge_weight_ratio = 0.05;
double transcript_min_expression = 15.0;
int min_hyper_edges_count = 20;
bool strand_reverse = false;

// for simulation
int simulation_num_vertices = 0;
int simulation_num_edges = 0;
int simulation_max_edge_weight = 0;

//// from command line
string algo;
string input_file;
string ref_file;
string output_file;

bool output_tex_files = false;
string fixed_gene_name = "";
int min_gtf_transcripts_num = 0;
bool fast_mode = true;

int print_parameters()
{
	// for bam file
	printf("parameters:\n");
	printf("min_flank_length = %d\n", min_flank_length);
	printf("min_bundle_gap = %d\n", min_bundle_gap);
	printf("min_subregion_gap = %d\n", min_subregion_gap);
	printf("min_num_hits_in_bundle = %d\n", min_num_hits_in_bundle);
	printf("min_splice_boundary_hits = %d\n", min_splice_boundary_hits);
	printf("min_mapping_quality = %d\n", min_mapping_quality);
	printf("min_average_overlap = %.2lf\n", min_average_overlap);
	printf("max_indel_ratio = %.2lf\n", max_indel_ratio);
	printf("min_subregion_length = %d\n", min_subregion_length);
	printf("min_max_region_overlap = %d\n", min_max_region_overlap);
	printf("min_region_coverage = %.2lf\n", min_region_coverage);
	printf("max_num_bundles = %d\n", max_num_bundles);
	printf("tail_coverage = %d\n", tail_coverage);
	printf("use_paired_end = %c\n", use_paired_end ? 'T' : 'F');
	printf("identify_slopes = %c\n", identify_slopes ? 'T' : 'F');
	printf("slope_bin_size = %d\n", slope_bin_size);
	printf("slope_min_bin_num = %d\n", slope_min_bin_num);
	printf("slope_min_score = %d\n", slope_min_score);
	printf("slope_extend_score = %d\n", slope_extend_score);
	printf("slope_flexible_bin_num = %d\n", slope_flexible_bin_num);
	printf("slope_acceptance_sigma = %.2lf\n", slope_acceptance_sigma);
	printf("average_slope_length = %d\n", average_slope_length);
	printf("average_read_length = %d\n", average_read_length);
	printf("pseudo_length_count = %d\n", pseudo_length_count);
	printf("strand_reverse = %c\n", strand_reverse ? 'T' : 'F');

	// for algorithm
	printf("join_min_reliability = %.2lf\n", join_min_reliability);
	printf("infer_min_reliability = %.2lf\n", infer_min_reliability);
	printf("infer_root_reliability = %.2lf\n", infer_root_reliability);
	printf("max_dp_table_size = %d\n", max_dp_table_size);
	printf("max_num_subsetsum_solutions = %d\n", max_num_subsetsum_solutions);
	printf("max_equation_error_ratio = %.2lf\n", max_equation_error_ratio);
	printf("max_gateway_error_ratio = %.2lf\n", max_gateway_error_ratio);
	printf("min_boundary_edge_weight_ratio = %.2lf\n", min_boundary_edge_weight_ratio);
	printf("transcript_min_expression = %.2lf\n", transcript_min_expression);
	printf("min_hyper_edges_count = %d\n", min_hyper_edges_count);
	printf("max_equations_each_iteration = %d\n", max_equations_each_iteration);

	// for simulation
	printf("simulation_num_vertices = %d\n", simulation_num_vertices);
	printf("simulation_num_edges = %d\n", simulation_num_edges);
	printf("simulation_max_edge_weight = %d\n", simulation_max_edge_weight);

	// for command
	printf("algo = %s\n", algo.c_str());
	printf("input_file = %s\n", input_file.c_str());
	printf("ref_file = %s\n", ref_file.c_str());
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
		else if(string(argv[i]) == "-g")
		{
			fixed_gene_name = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-t")
		{
			output_tex_files = true;
		}
		else if(string(argv[i]) == "-p")
		{
			use_paired_end = true;
		}
		else if(string(argv[i]) == "-RF")
		{
			strand_reverse = true;
		}
		else if(string(argv[i]) == "-x")
		{
			transcript_min_expression = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-R")
		{
			// default setting for real dataset
			min_splice_boundary_hits = 2;
			min_num_hits_in_bundle = 20;
			min_bundle_gap = 50;
			min_mapping_quality = 1;
		}
	}

	return b;
}

