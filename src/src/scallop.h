#ifndef __SCALLOP_H__
#define __SCALLOP_H__

#include "sgraph.h"

class scallop
{
public:
	vector<sgraph> sgraphs;

public:
	scallop();
	~scallop();

public:
	int process(const char *bam_file);
	int load(const char *bam_file);
	int solve();
};

#endif
