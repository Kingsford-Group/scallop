/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "config.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>

using namespace std;

config::config(){
	//// parameters
	// for bam file and reads
	min_flank_length = 3;
	max_edit_distance = 10;
	min_bundle_gap = 50;
	min_num_hits_in_bundle = 20;
	min_mapping_quality = 1;
	min_splice_boundary_hits = 1;
	use_second_alignment = false;
	uniquely_mapped_only = false;
	library_type = EMPTY;

	// for preview
	max_preview_reads = 2000000;
	max_preview_spliced_reads = 50000;
	min_preview_spliced_reads = 10000;
	preview_infer_ratio = 0.95;
	preview_only = false;

	// for identifying subgraphs
	min_subregion_gap = 3;
	min_subregion_overlap = 1.5;
	min_subregion_length = 15;

	// for revising/decomposing splice graph
	max_intron_contamination_coverage = 2.0;
	min_surviving_edge_weight = 1.5;
	max_decompose_error_ratio[0] = 0.33;
	max_decompose_error_ratio[1] = 0.05;
	max_decompose_error_ratio[2] = 0.0;
	max_decompose_error_ratio[3] = 0.25;
	max_decompose_error_ratio[4] = 0.30;
	max_decompose_error_ratio[5] = 0.0;
	max_decompose_error_ratio[6] = 1.1;

	// for selecting paths
	min_transcript_coverage = 1.01;
	min_transcript_coverage_ratio = 0.005;
	min_single_exon_coverage = 20;
	min_transcript_numreads = 20;
	min_transcript_length_base = 150;
	min_transcript_length_increase = 50;
	min_exon_length = 20;
	max_num_exons = 1000;

	// for subsetsum and router
	max_dp_table_size = 10000;
	min_router_count = 1;

	// for simulation
	simulation_num_vertices = 0;
	simulation_num_edges = 0;
	simulation_max_edge_weight = 0;

	// input and output
	algo = "scallop";
	input_file;
	ref_file;
	ref_file1;
	ref_file2;
	output_file;

	// for controling
	output_tex_files = false;
	fixed_gene_name = "";
	batch_bundle_size = 100;
	verbose = 1;
	version = "v0.10.2b";
}


void config::update_from_file(char * fname)
{
	ifstream f(fname);
	if(f){
		string arg = "";
		while(f.peek() != EOF){
			f >> arg;
			if(arg == "min_flank_length"){
				f >> min_flank_length;
			}else if(arg == "max_edit_distance"){
				f >> max_edit_distance;
			}else if(arg == "min_bundle_gap"){
				f >> min_bundle_gap;
			}else if(arg == "min_num_hits_in_bundle"){
				f >> min_num_hits_in_bundle;
			}else if(arg == "min_mapping_quality"){
				f >> min_mapping_quality;
			}else if(arg == "min_splice_boundary_hits"){
				f >> min_splice_boundary_hits;
			}else if(arg == "min_subregion_gap"){
				f >> min_subregion_gap;
			}else if(arg == "min_subregion_length"){
				f >> min_subregion_length;
			}else if(arg == "min_subregion_overlap"){
				f >> min_subregion_overlap;
			}else if(arg == "min_surviving_edge_weight"){
				f >> min_surviving_edge_weight;
			}else if(arg == "max_intron_contamination_coverage"){
				f >> max_intron_contamination_coverage;
			}else if(arg == "min_transcript_coverage"){
				f >> min_transcript_coverage;
				if(fabs(min_transcript_coverage - 1.0) < 0.01) min_transcript_coverage = 1.01;
			}else if(arg == "min_transcript_coverage_ratio"){
				f >> min_transcript_coverage_ratio;
			}else if(arg == "min_single_exon_coverage"){
				f >> min_single_exon_coverage;
			}else if(arg == "min_transcript_numreads"){
				f >> min_transcript_numreads;
			}else if(arg == "min_transcript_length_base"){
				f >> min_transcript_length_base;
			}else if(arg == "min_transcript_length_increase"){
				f >> min_transcript_length_increase;
			}else if(arg == "min_exon_length"){
				f >> min_exon_length;
			}else if(arg == "max_num_exons"){
				f >> max_num_exons;
			}else if(arg == "max_dp_table_size"){
				f >> max_dp_table_size;
			}else if(arg == "min_router_count"){
				f >> min_router_count;
			}else if(arg == "use_second_alignment"){
				string s;
				f >> s;
				std::transform(s.begin(), s.end(), s.begin(), ::tolower);

				if(s == "true" or s == "1") use_second_alignment = true;
				else if(s == "false" or s == "0") use_second_alignment = false;
				else{
					std::cerr << "'" << s << "' not recognzed as a boolean for use_second_alignment" << std::endl;
					exit(64);
				}
			}else if(arg == "uniquely_mapped_only"){
				string s;
				f >> s;
				std::transform(s.begin(), s.end(), s.begin(), ::tolower);

				if(s == "true" or s == "1") uniquely_mapped_only = true;
				else if(s == "false" or s == "0") uniquely_mapped_only = false;
				else{
					std::cerr << "'" << s << "' not recognzed as a boolean for uniquely_mapped_only" << std::endl;
					exit(64);
				}
			}else if(arg == "batch_bundle_size"){
				f >> batch_bundle_size;
			}else{
				cerr << "Unknown option in configuration file" << endl;
				cerr << "File: " << fname << endl;
				cerr << "Argument: " << arg << endl;
				exit(64);
			}
		}
	}else{
		cerr << "Cannot open configuration file " << fname << endl;
		cerr << "Error code: " << strerror(errno);
		exit(126);
	}
}

int config::parse_arguments(int argc, const char ** argv)
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
		else if(string(argv[i]) == "--version")
		{
			printf("%s\n", version.c_str());
			exit(0);
		}
		else if(string(argv[i]) == "--help")
		{
			print_copyright();
			print_help();
			printf("\n");
			print_logo();
			exit(0);
		}
		else if(string(argv[i]) == "--min_flank_length")
		{
			min_flank_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_edit_distance")
		{
			max_edit_distance = atoi(argv[i + 1]);
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
		else if(string(argv[i]) == "--max_preview_spliced_reads")
		{
			max_preview_spliced_reads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_preview_spliced_reads")
		{
			min_preview_spliced_reads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--preview")
		{
			preview_only = true;
		}
		else if(string(argv[i]) == "--max_preview_reads")
		{
			max_preview_reads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--preview_infer_ratio")
		{
			preview_infer_ratio = atof(argv[i + 1]);
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
		else if(string(argv[i]) == "--min_surviving_edge_weight")
		{
			min_surviving_edge_weight = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_intron_contamination_coverage")
		{
			max_intron_contamination_coverage = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_coverage")
		{
			min_transcript_coverage = atof(argv[i + 1]);
			i++;
			if(fabs(min_transcript_coverage - 1.0) < 0.01) min_transcript_coverage = 1.01;
		}
		else if(string(argv[i]) == "--min_transcript_coverage_ratio")
		{
			min_transcript_coverage_ratio = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_single_exon_coverage")
		{
			min_single_exon_coverage = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_numreads")
		{
			min_transcript_numreads = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_length_base")
		{
			min_transcript_length_base = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_length_increase")
		{
			min_transcript_length_increase = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_exon_length")
		{
			min_exon_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_num_exons")
		{
			max_num_exons = atoi(argv[i + 1]);
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
			if(s == "empty") library_type = EMPTY;
			if(s == "unstranded") library_type = UNSTRANDED;
			if(s == "first") library_type = FR_FIRST;
			if(s == "second") library_type = FR_SECOND;
			i++;
		}
		else if(string(argv[i]) == "--use_second_alignment")
		{
			string s(argv[i + 1]);
			std::transform(s.begin(), s.end(), s.begin(), ::tolower);

			if(s == "true" or s == "1") use_second_alignment = true;
			else if(s == "false" or s == "0") use_second_alignment = false;
			else{
				std::cerr << "'" << s << "' not recognzed as a boolean for use_second_alignment" << std::endl;
				exit(64);
			}
			i++;
		}
		else if(string(argv[i]) == "--uniquely_mapped_only")
		{
			string s(argv[i + 1]);
			std::transform(s.begin(), s.end(), s.begin(), ::tolower);

			if(s == "true" or s == "1") uniquely_mapped_only = true;
			else if(s == "false" or s == "0") uniquely_mapped_only = false;
			else{
				std::cerr << "'" << s << "' not recognzed as a boolean for uniquely_mapped_only" << std::endl;
				exit(64);
			}
			i++;
		}
		else if(string(argv[i]) == "--verbose")
		{
			verbose = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--batch_bundle_size")
		{
			batch_bundle_size = atoi(argv[i + 1]);
			i++;
		}
	}

	if(min_surviving_edge_weight < 0.1 + min_transcript_coverage)
	{
		min_surviving_edge_weight = 0.1 + min_transcript_coverage;
	}

	// verify arguments
	if(input_file == "")
	{
		printf("error: input-file is missing.\n");
		exit(0);
	}

	if(output_file == "" && preview_only == false)
	{
		printf("error: output-file is missing.\n");
		exit(0);
	}

	return 0;
}

int config::print_parameters()
{
	printf("parameters:\n");

	// for bam file and reads
	printf("min_flank_length = %d\n", min_flank_length);
	printf("max_edit_distance = %d\n", max_edit_distance);
	printf("min_bundle_gap = %d\n", min_bundle_gap);
	printf("min_num_hits_in_bundle = %d\n", min_num_hits_in_bundle);
	printf("min_mapping_quality = %d\n", min_mapping_quality);
	printf("min_splice_boundary_hits = %d\n", min_splice_boundary_hits);

	// for preview
	printf("preview_only = %c\n", preview_only ? 'T' : 'F');
	printf("max_preview_reads = %d\n", max_preview_reads);
	printf("max_preview_spliced_reads = %d\n", max_preview_spliced_reads);
	printf("min_preview_spliced_reads = %d\n", min_preview_spliced_reads);
	printf("preview_infer_ratio = %.3lf\n", preview_infer_ratio);

	// for identifying subgraphs
	printf("min_subregion_gap = %d\n", min_subregion_gap);
	printf("min_subregion_length = %d\n", min_subregion_length);
	printf("min_subregion_overlap = %.2lf\n", min_subregion_overlap);

	// for splice graph
	printf("max_intron_contamination_coverage = %.2lf\n", max_intron_contamination_coverage);
	printf("min_surviving_edge_weight = %.2lf\n", min_surviving_edge_weight);
	printf("min_transcript_coverage = %.2lf\n", min_transcript_coverage);
	printf("min_transcript_coverage_ratio = %.2lf\n", min_transcript_coverage_ratio);
	printf("min_single_exon_coverage = %.2lf\n", min_single_exon_coverage);
	printf("min_transcript_numreads = %.2lf\n", min_transcript_numreads);
	printf("min_transcript_length_base = %d\n", min_transcript_length_base);
	printf("min_transcript_length_increase = %d\n", min_transcript_length_increase);
	printf("max_num_exons = %d\n", max_num_exons);

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
	printf("library_type = %d\n", library_type);
	printf("output_tex_files = %c\n", output_tex_files ? 'T' : 'F');
	printf("fixed_gene_name = %s\n", fixed_gene_name.c_str());
	printf("use_second_alignment = %c\n", use_second_alignment ? 'T' : 'F');
	printf("uniquely_mapped_only = %c\n", uniquely_mapped_only ? 'T' : 'F');
	printf("verbose = %d\n", verbose);
	printf("batch_bundle_size = %d\n", batch_bundle_size);

	printf("\n");

	return 0;
}

int config::print_command_line(int argc, const char ** argv)
{
	printf("command line: ");
	for(int i = 0; i < argc; i++)
	{
		printf("%s ", argv[i]);
	}
	printf("\n");
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
	printf("\n");
	printf("Usage: scallop -i <bam-file> -o <gtf-file> [options]\n");
	printf("\n");
	printf("Options:\n");
	printf(" %-42s  %s\n", "--help",  "print usage of Scallop and exit");
	printf(" %-42s  %s\n", "--version",  "print current version of Scallop and exit");
	printf(" %-42s  %s\n", "--verbose <0, 1, 2>",  "0: quiet; 1: one line for each graph; 2: with details, default: 1");
	printf(" %-42s  %s\n", "--library_type <first, second, unstranded>",  "library type of the sample, default: unstranded");
	printf(" %-42s  %s\n", "--min_transcript_coverage <float>",  "minimum coverage required for a multi-exon transcript, default: 1.01");
	printf(" %-42s  %s\n", "--min_single_exon_coverage <float>",  "minimum coverage required for a single-exon transcript, default: 20");
	printf(" %-42s  %s\n", "--min_transcript_length_increase <integer>",  "default: 50");
	printf(" %-42s  %s\n", "--min_transcript_length_base <integer>",  "default: 250, minimum length of a transcript would be");
	printf(" %-42s  %s\n", "",  "--min_transcript_length_base + --min_transcript_length_increase * num-of-exons");
	printf(" %-42s  %s\n", "--min_mapping_quality <integer>",  "ignore reads with mapping quality less than this value, default: 1");
	printf(" %-42s  %s\n", "--min_bundle_gap <integer>",  "minimum distances required to start a new bundle, default: 50");
	printf(" %-42s  %s\n", "--min_num_hits_in_bundle <integer>",  "minimum number of reads required in a bundle, default: 20");
	printf(" %-42s  %s\n", "--min_flank_length <integer>",  "minimum match length in each side for a spliced read, default: 3");
	printf(" %-42s  %s\n", "--min_splice_bundary_hits <integer>",  "minimum number of spliced reads required for a junction, default: 1");
	return 0;
}

int config::print_copyright()
{
	printf("Scallop %s (c) 2017 Mingfu Shao, Carl Kingsford, and Carnegie Mellon University\n", version.c_str());
	return 0;
}
