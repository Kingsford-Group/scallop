#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "util.h"
#include <stdint.h>
#include <map>
#include <sstream>

using namespace std;

// constants
#define LEFT_SPLICE 1
#define RIGHT_SPLICE 2
#define LEFT_BOUNDARY 3
#define RIGHT_BOUNDARY 4
#define START_BOUNDARY 5
#define END_BOUNDARY 6

#define TRIVIAL 0
#define NORMAL 1

#define SLOPE5END 0
#define SLOPE3END 1

#define SLOPE_MARGIN 0
#define SLOPE_MIDDLE 1

// pre-defined parameters
#define MAX_NUM_CIGAR 7
#define MIN_LEN_FLANK 1

// user-defined parameters
extern int32_t min_bundle_gap;
extern int min_num_hits_in_bundle;
extern int32_t min_splice_boundary_hits;
extern uint32_t min_max_splice_boundary_qual;
extern int32_t average_read_length;
extern int32_t average_slope_length;
extern double min_average_overlap;
extern int min_max_region_overlap;
extern double min_region_coverage;
extern int max_num_bundles;
extern int max_dp_table_size;
extern int max_num_subsetsum_solutions;
extern double max_equation_error_ratio;
extern int slope_bin_size;
extern int slope_min_bin_num;
extern int slope_std_bin_num;
extern int slope_min_distance;
extern int min_slope_score;
extern int pseudo_length_count;
extern double min_boundary_edge_weight_ratio;

extern string algo;
extern string input_file;
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

// load parameters
int load_config(const char * conf_file);

#endif
