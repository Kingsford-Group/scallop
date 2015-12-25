#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <stdint.h>

// constants
#define SPLICE_BOUNDARY 0
#define START_BOUNDARY 1
#define END_BOUNDARY 2

// pre-defined parameters
#define MAX_NUM_CIGAR 7
#define MIN_LEN_FLANK 1

// user-defined parameters
extern int32_t min_bundle_gap;

extern int32_t min_splice_boundary_hits;
extern int32_t min_start_boundary_hits;
extern int32_t min_end_boundary_hits;
extern uint32_t min_max_splice_boundary_qual;
extern uint32_t min_max_start_boundary_qual;
extern uint32_t min_max_end_boundary_qual;
extern int32_t hits_window_size;

// load parameters
int load_config(const char * conf_file);

#endif
