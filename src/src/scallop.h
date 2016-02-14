#ifndef __SCALLOP_H__
#define __SCALLOP_H__

#include "bundle.h"

class scallop
{
public:
	vector<bundle> bundles;

public:
	scallop();
	~scallop();

public:
	int process(const string &file);
	int load(const char *bam_file);
};

#endif
