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
bool use_second_alignment = false;
int library_type = UNSTRAND;

// for identifying subgraphs
double max_indel_ratio = 0.2;
int32_t min_subregion_gap = 3;
double min_subregion_overlap = 1.5;
int32_t min_subregion_length = 15;

// for revising/decomposing splice graph
double min_splice_edge_weight = 1.5;
double max_small_error_ratio = 0.33;
double max_split_error_ratio = 0.25;
double max_decompose_error_ratio = 0.01;

// for selecting paths
double min_transcript_coverage = 0.9;
double min_transcript_numreads = 20;
int min_transcript_length = 200;

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
int32_t average_read_length = 100;
bool output_tex_files = false;
string fixed_gene_name = "";

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
	printf("min_splice_edge_weight = %.2lf\n", min_splice_edge_weight);
	printf("min_transcript_coverage = %.2lf\n", min_transcript_coverage);
	printf("min_transcript_numreads = %.2lf\n", min_transcript_numreads);
	printf("min_transcript_length = %d\n", min_transcript_length);
	printf("max_small_error_ratio = %.2lf\n", max_small_error_ratio);
	printf("max_split_error_ratio = %.2lf\n", max_split_error_ratio);
	printf("max_decompose_error_ratio = %.2lf\n", max_decompose_error_ratio);

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
	printf("average_read_length = %d\n", average_read_length);
	printf("library_type = %d\n", library_type);
	printf("output_tex_files = %c\n", output_tex_files ? 'T' : 'F');
	printf("fixed_gene_name = %s\n", fixed_gene_name.c_str());
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
		else if(string(argv[i]) == "--min_splice_edge_weight")
		{
			min_splice_edge_weight = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_coverage")
		{
			min_transcript_coverage = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_numreads")
		{
			min_transcript_numreads = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_length")
		{
			min_transcript_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_split_error_ratio")
		{
			max_split_error_ratio = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_small_error_ratio")
		{
			max_small_error_ratio = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_decompose_error_ratio")
		{
			max_decompose_error_ratio = atof(argv[i + 1]);
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
		else if(string(argv[i]) == "--library_type")
		{
			string s(argv[i + 1]);
			if(s == "unstrand") library_type = UNSTRAND;
			if(s == "first") library_type = FR_FIRST;
			if(s == "second") library_type = FR_SECOND;
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

	// TODO
	//min_splice_edge_weight = min_transcript_coverage;

	return 0;
}

int print_logo()
{
	printf("      ___           ___           ___                                       ___           ___    \n");
	printf("     /  /\\         /  /\\         /  /\\                                     /  /\\         /  /\\   \n");
	printf("    /  /:/_       /  /:/        /  /::\\                                   /  /::\\       /  /::\\  \n");
	printf("   /  /:/ /\\     /  /:/        /  /:/\\:\\    ___     ___   ___     ___    /  /:/\\:\\     /  /:/\\:\\ \n");
	printf("  /  /:/ /::\\   /  /:/  ___   /  /:/~/::\\  /__/\\   /  /\\ /__/\\   /  /\\  /  /:/  \\:\\   /  /:/~/:/ \n");
	printf(" /__/:/ /:/\\:\\ /__/:/  /  /\\ /__/:/ /:/\\:\\ \\  \\:\\ /  /:/ \\  \\:\\ /  /:/ /__/:/ \\__\\:\\ /__/:/ /:/  \n");
	printf(" \\  \\:\\/:/~/:/ \\  \\:\\ /  /:/ \\  \\:\\/:/__\\/  \\  \\:\\  /:/   \\  \\:\\  /:/  \\  \\:\\ /  /:/ \\  \\:\\/:/   \n");
	printf("  \\  \\::/ /:/   \\  \\:\\  /:/   \\  \\::/        \\  \\:\\/:/     \\  \\:\\/:/    \\  \\:\\  /:/   \\  \\::/    \n");
	printf("   \\__\\/ /:/     \\  \\:\\/:/     \\  \\:\\         \\  \\::/       \\  \\::/      \\  \\:\\/:/     \\  \\:\\    \n");
	printf("     /__/:/       \\  \\::/       \\  \\:\\         \\__\\/         \\__\\/        \\  \\::/       \\  \\:\\   \n");
	printf("     \\__\\/         \\__\\/         \\__\\/                                     \\__\\/         \\__\\/   \n");
	printf("\n");

	return 0;
}

int print_help()
{
	printf("typical command line: ./scallop -i bam-file -o gtf-file\n");
	return 0;
}
