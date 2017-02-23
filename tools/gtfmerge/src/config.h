#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <stdint.h>
#include <map>
#include <sstream>

using namespace std;

extern int algo;
extern double min_transcript_length;
extern bool multiple_exon;

int parse_parameters(int argc, const char ** argv);

#endif
