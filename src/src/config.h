/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "util.h"
#include <stdint.h>
#include <map>
#include <sstream>

using namespace std;

//// constants
#define MAX_NUM_CIGAR 7

#define START_BOUNDARY 1
#define END_BOUNDARY 2
#define LEFT_SPLICE 3
#define RIGHT_SPLICE 4
#define LEFT_RIGHT_SPLICE 5
#define MIDDLE_CUT 6

#define TRIVIAL 0
#define NORMAL 1

// five types for decomposition
#define SMALLEST_EDGE 0
#define NEGLIGIBLE_EDGE 1
#define SPLITTABLE_SIMPLE 2
#define SPLITTABLE_HYPER 3
#define UNSPLITTABLE_SINGLE 4
#define UNSPLITTABLE_MULTIPLE 5
#define TRIVIAL_VERTEX 6

#define UNSTRANDED 0
#define FR_FIRST 1
#define FR_SECOND 2

//// parameters
// for bam file and reads
extern int min_flank_length;
extern int max_edit_distance;
extern int32_t min_bundle_gap;
extern int min_num_hits_in_bundle;
extern uint32_t min_mapping_quality;
extern int32_t min_splice_boundary_hits;
extern bool uniquely_mapped_only;
extern bool use_second_alignment;

// for identifying subgraphs
extern int32_t min_subregion_gap;
extern double min_subregion_overlap;
extern int32_t min_subregion_length;
extern int min_subregion_ladders;

// for subsetsum and router
extern int max_dp_table_size;
extern int min_router_count;

// for splice graph
extern double max_intron_contamination_coverage;
extern double min_surviving_edge_weight;
extern double max_decompose_error_ratio[7];
extern double min_transcript_numreads;
extern double min_transcript_coverage;
extern double min_single_exon_coverage;
extern double min_transcript_coverage_ratio; 
extern int min_transcript_length_base;
extern int min_transcript_length_increase;
extern int min_exon_length;
extern int max_num_exons;

// for simulation
extern int simulation_num_vertices;
extern int simulation_num_edges;
extern int simulation_max_edge_weight;

// input and output
extern string algo;
extern string input_file;
extern string ref_file;
extern string ref_file1;
extern string ref_file2;
extern string output_file;

// for controling
extern bool output_tex_files;
extern string fixed_gene_name;
extern int max_num_bundles;
extern int library_type;
extern int min_gtf_transcripts_num;
extern int batch_bundle_size;
extern int verbose;
extern string version;

// parse arguments
int print_command_line(int argc, const char ** argv);
int parse_arguments(int argc, const char ** argv);
int print_parameters();
int print_copyright();
int print_logo();
int print_help();

#endif
