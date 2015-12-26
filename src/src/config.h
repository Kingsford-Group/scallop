#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <stdint.h>

// constants
#define SPLICE_BOUNDARY 1
#define LEFT_BOUNDARY 2
#define RIGHT_BOUNDARY 3

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
extern int32_t hits_window_size;
extern double min_boundary_score;

// load parameters
int load_config(const char * conf_file);

#endif
