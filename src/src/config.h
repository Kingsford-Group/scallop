#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <stdint.h>

// pre-defined parameters
#define MAX_NUM_CIGAR 7
#define MIN_LEN_FLANK 1

// user-defined parameters
static int32_t min_bundle_gap;

// load parameters
int load_config(const char * conf_file);

#endif
