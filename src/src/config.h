#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <stdint.h>
#include <map>

using namespace std;

// macros: using int64_t for two int32_t
#define pack(x, y) (int64_t)((((int64_t)(x)) << 32) | ((int64_t)(y)))
#define high32(x) (int32_t)((x) >> 32)
#define low32(x) (int32_t)(((x) << 32) >> 32)


// definitions
typedef map<int32_t, int> MPI;
typedef pair<int32_t, int> PPI;


// constants
#define LEFT_SPLICE 1
#define RIGHT_SPLICE 2
#define LEFT_BOUNDARY 3
#define RIGHT_BOUNDARY 4
#define START_BOUNDARY 5
#define END_BOUNDARY 6

// pre-defined parameters
#define MAX_NUM_CIGAR 7
#define MIN_LEN_FLANK 1

// user-defined parameters
extern int32_t min_bundle_gap;

extern int32_t min_splice_boundary_hits;
extern int32_t min_left_boundary_hits;
extern int32_t min_right_boundary_hits;
extern uint32_t min_max_splice_boundary_qual;
extern uint32_t min_max_left_boundary_qual;
extern uint32_t min_max_right_boundary_qual;
extern int32_t average_read_length;
extern uint32_t min_boundary_score;
extern int32_t ascending_step;
extern int32_t descending_step;
extern uint32_t min_ascending_score;
extern uint32_t min_descending_score;
extern int num_sample_positions;
extern double min_average_overlap;
extern int min_max_region_overlap;
extern double min_region_coverage;

// load parameters
int load_config(const char * conf_file);

#endif
