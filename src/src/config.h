#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "util.h"
#include <stdint.h>
#include <map>
#include <sstream>

using namespace std;

//// constants
#define START_BOUNDARY 1
#define END_BOUNDARY 2
#define LEFT_SPLICE 3
#define RIGHT_SPLICE 4
#define LEFT_RIGHT_SPLICE 5
#define MIDDLE_CUT 6

#define TRIVIAL 0
#define NORMAL 1

#define MAX_NUM_CIGAR 7

#define SPLITABLE 1
#define SINGLE 2
#define MULTIPLE 3

//// parameters
// for bam file and reads
extern int min_flank_length;
extern int32_t min_bundle_gap;
extern int min_num_hits_in_bundle;
extern uint32_t min_mapping_quality;
extern int32_t min_splice_boundary_hits;

// for identifying subgraphs
extern double max_indel_ratio;
extern int32_t min_subregion_gap;
extern double min_subregion_overlap;
extern int32_t min_subregion_length;
extern int min_subregion_ladders;

// for subsetsum and router
extern int max_dp_table_size;
extern int min_router_count;

// for splice graph
extern double min_inner_vertex_weight;
extern double min_inner_boundary_weight;
extern double min_spanning_edge_weight;
extern double min_splice_edge_weight;
extern double max_split_error_ratio;
extern double min_transcript_coverage;
extern double min_splice_graph_coverage;
extern int min_junction_count;
extern double smallest_edge_ratio_scalor1;
extern double smallest_edge_ratio_scalor2;
extern double max_removable_weight;
extern double min_removable_weight;
extern bool extend_isolated_boundary;

// for identifying new boundaries
extern bool identify_extra_boundary;
extern int min_boundary_score;
extern int min_boundary_length;
extern double min_boundary_ave_ratio;
extern double min_boundary_sigma;

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
extern int32_t average_read_length;
extern bool output_tex_files;
extern string fixed_gene_name;
extern int max_num_bundles;
extern bool strand_reverse;
extern int min_gtf_transcripts_num;
extern bool fast_mode;
extern bool use_second_alignment;

// parse arguments
int print_command_line(int argc, const char ** argv);
int parse_arguments(int argc, const char ** argv);
int print_parameters();

#endif
