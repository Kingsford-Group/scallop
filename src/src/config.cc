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
int min_subregion_ladders = 0;

// for splice graph
double min_inner_vertex_weight = 10;
double min_inner_boundary_weight = 4.0;
double min_splice_edge_weight = 3.5;
double min_spanning_edge_weight = 5.5;
double max_split_error_ratio = 0.15;
double max_decompose_error_ratio = 0.05;
double min_transcript_coverage = 10.0;
double min_splice_graph_coverage = 20.0;
int min_junction_count = 5;
double smallest_edge_ratio_scalor1 = 0.2;
double smallest_edge_ratio_scalor2 = 1.0;
double min_removable_weight = 5.0;
double max_removable_weight = 50.0;
bool extend_isolated_boundary = true;

// for identifying new boundaries
bool identify_extra_boundary = false;
int min_boundary_length = 80;
int min_boundary_score = 100;
double min_boundary_ave_ratio = 3.0;
double min_boundary_sigma = 5.0;

// for subsetsum and router
int max_dp_table_size = 10000;
int min_router_count = 1;

// for simulation
int simulation_num_vertices = 0;
int simulation_num_edges = 0;
int simulation_max_edge_weight = 0;

// input and output
string algo = "scallop";
string input_file;
string ref_file;
string ref_file1;
string ref_file2;
string output_file;

// for controling
int max_num_bundles = -1;
int32_t average_read_length = 100;
bool strand_reverse = false;
bool output_tex_files = false;
string fixed_gene_name = "";
int min_gtf_transcripts_num = 0;
bool fast_mode = true;
bool use_second_alignment = false;

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
	printf("min_subregion_ladders = %d\n", min_subregion_ladders);

	// for splice graph
	printf("min_splice_edge_weight = %.2lf\n", min_splice_edge_weight);
	printf("min_spanning_edge_weight = %.2lf\n", min_spanning_edge_weight);
	printf("min_inner_vertex_weight = %.2lf\n", min_inner_vertex_weight);
	printf("min_inner_boundary_weight = %.2lf\n", min_inner_boundary_weight);
	printf("min_transcript_coverage = %.2lf\n", min_transcript_coverage);
	printf("min_splice_graph_coverage = %.2lf\n", min_splice_graph_coverage);
	printf("min_junction_count = %d\n", min_junction_count);
	printf("max_split_error_ratio = %.2lf\n", max_split_error_ratio);
	printf("max_decompose_error_ratio = %.2lf\n", max_decompose_error_ratio);
	printf("smallest_edge_ratio_scalor1 = %.2lf\n", smallest_edge_ratio_scalor1);
	printf("smallest_edge_ratio_scalor2 = %.2lf\n", smallest_edge_ratio_scalor2);
	printf("min_removable_weight = %.2lf\n", min_removable_weight);
	printf("max_removable_weight = %.2lf\n", max_removable_weight);
	printf("extend_isolated_bounary = %c\n", extend_isolated_boundary ? 'T' : 'F');

	// for identifying new boundaries
	printf("identify_extra_boundary = %c\n", identify_extra_boundary ? 'T' : 'F');
	printf("min_boundary_length = %d\n", min_boundary_length);
	printf("min_boundary_score = %d\n", min_boundary_score);
	printf("min_boundary_ave_ratio = %.2lf\n", min_boundary_ave_ratio);
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
	printf("output_tex_files = %c\n", output_tex_files ? 'T' : 'F');
	printf("fixed_gene_name = %s\n", fixed_gene_name.c_str());
	printf("min_gtf_transcripts_num = %d\n", min_gtf_transcripts_num);
	printf("fast_mode = %c\n", fast_mode ? 'T' : 'F');
	printf("use_second_alignment = %c\n", use_second_alignment ? 'T' : 'F');

	printf("\n");

	return 0;
}

int print_command_line(int argc, const char ** argv)
{
	printf("command line: ");
	for(int i = 0; i < argc; i++)
	{
		printf("%s ", argv[i]);
	}
	printf("\n");
	return 0;
}

int parse_arguments(int argc, const char ** argv)
{
	for(int i = 1; i < argc; i++)
	{
		// necessary ones
		if(string(argv[i]) == "-i")
		{
			input_file = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-o")
		{
			output_file = string(argv[i + 1]);
			i++;
		}

		// internal use
		else if(string(argv[i]) == "-a")
		{
			algo = string(argv[i + 1]);
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

		// user specified
		else if(string(argv[i]) == "--min_flank_length")
		{
			min_flank_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_bundle_gap")
		{
			min_bundle_gap = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_num_hits_in_bundle")
		{
			min_num_hits_in_bundle = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_mapping_quality")
		{
			min_mapping_quality = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_splice_boundary_hits")
		{
			min_splice_boundary_hits = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_indel_ratio")
		{
			max_indel_ratio = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_gap")
		{
			min_subregion_gap = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_length")
		{
			min_subregion_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_overlap")
		{
			min_subregion_overlap = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--identify_extra_boundary")
		{
			string s(argv[i + 1]);
			if(s == "true") identify_extra_boundary = true;
			else identify_extra_boundary = false;
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_ladders")
		{
			min_subregion_ladders = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_splice_edge_weight")
		{
			min_splice_edge_weight = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_spanning_edge_weight")
		{
			min_spanning_edge_weight = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_inner_vertex_weight")
		{
			min_inner_vertex_weight = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_inner_boundary_weight")
		{
			min_inner_boundary_weight = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_coverage")
		{
			min_transcript_coverage = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_splice_graph_coverage")
		{
			min_splice_graph_coverage = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_junction_count")
		{
			min_junction_count = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_split_error_ratio")
		{
			max_split_error_ratio = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_decompose_error_ratio")
		{
			max_decompose_error_ratio = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--smallest_edge_ratio_scalor1")
		{
			smallest_edge_ratio_scalor1 = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--smallest_edge_ratio_scalor2")
		{
			smallest_edge_ratio_scalor2 = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_removable_weight")
		{
			min_removable_weight = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_removable_weight")
		{
			max_removable_weight = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--extend_isolated_boundary")
		{
			string s(argv[i + 1]);
			if(s == "true") extend_isolated_boundary = true;
			else extend_isolated_boundary = false;
			i++;
		}
		else if(string(argv[i]) == "--min_boundary_length")
		{
			min_boundary_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_boundary_score")
		{
			min_boundary_score = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_boundary_ave_ratio")
		{
			min_boundary_ave_ratio = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_boundary_sigma")
		{
			min_boundary_sigma = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_dp_table_size")
		{
			max_dp_table_size = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_router_count")
		{
			min_router_count = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--strand_reverse")
		{
			string s(argv[i + 1]);
			if(s == "true") strand_reverse = true;
			else strand_reverse = false;
			i++;
		}
		else if(string(argv[i]) == "--use_second_alignment")
		{
			string s(argv[i + 1]);
			if(s == "true") use_second_alignment = true;
			else use_second_alignment = false;
			i++;
		}
	}

	return 0;
}
