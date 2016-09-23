#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "util.h"
#include <stdint.h>
#include <map>
#include <sstream>

using namespace std;

// constants
#define START_BOUNDARY 1
#define END_BOUNDARY 2
#define LEFT_SPLICE 3
#define RIGHT_SPLICE 4
#define LEFT_RIGHT_SPLICE 5
#define MIDDLE_CUT 6

#define TRIVIAL 0
#define NORMAL 1

#define SLOPE5END START_BOUNDARY
#define SLOPE3END END_BOUNDARY

#define SLOPE_MARGIN 0
#define SLOPE_MIDDLE 1

#define MAX_NUM_CIGAR 7

// user-defined parameters
extern int min_flank_length;
extern int min_clip_length;
extern int32_t min_bundle_gap;
extern int32_t min_subregion_gap;
extern double min_subregion_overlap;
extern double min_average_overlap;
extern int min_num_hits_in_bundle;
extern int32_t min_splice_boundary_hits;
extern uint32_t min_mapping_quality;
extern int32_t average_read_length;
extern double max_indel_ratio;
extern int32_t min_subregion_length;
extern int min_max_region_overlap;
extern double min_region_coverage;
extern int max_num_bundles;
extern int max_dp_table_size;
extern int max_num_subsetsum_solutions;
extern int max_equations_each_iteration;
extern double max_equation_error_ratio;
extern double max_split_error_ratio;
extern double max_router_error_ratio;
extern int min_router_count;
extern int tail_coverage;
extern int slope_bin_size;
extern int slope_bin_num;
extern int slope_min_score;
extern double slope_min_sigma;
extern bool use_paired_end;
extern bool ignore_single_exon_transcripts;
extern double max_ignorable_edge_weight;
extern double max_removable_edge_weight;
extern double min_edge_weight;
extern double min_transcript_coverage;
extern double min_splice_graph_coverage;
extern int min_boundary_score;
extern int min_boundary_length;
extern double min_boundary_sigma;
extern int32_t partial_exon_length;

extern double join_min_reliability;
extern double infer_min_reliability;
extern double infer_root_reliability;
extern int pseudo_length_count;
extern double min_boundary_edge_weight_ratio;
extern double transcript_min_expression;
extern bool strand_reverse;

extern int min_hyper_edges_count;

extern string algo;
extern string input_file;
extern string ref_file;
extern string ref_file1;
extern string ref_file2;
extern string output_file;
extern bool output_tex_files;
extern string fixed_gene_name;

extern int min_gtf_transcripts_num;
extern int simulation_num_vertices;
extern int simulation_num_edges;
extern int simulation_max_edge_weight;
extern bool fast_mode;

// parse arguments
bool parse_arguments(int argc, const char ** argv);
int print_parameters();

// load parameters
int load_config(const char * conf_file);

#endif
