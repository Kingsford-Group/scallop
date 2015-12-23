#ifndef __SCALLOP_H__
#define __SCALLOP_H__

#include "config.h"
#include "bundle.h"

class scallop
{
public:
	config *conf;
	vector<bundle> bundles;

public:
	scallop(config *_conf);
	~scallop();

public:
	int process(const char * bam_file);
};

#endif
