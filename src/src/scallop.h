#ifndef __SCALLOP_H__
#define __SCALLOP_H__

#include "config.h"
#include "bundle.h"

class scallop
{
public:
	config *conf;
	bundle **bundles;
	int n_bundles;
	int m_bundles;

public:
	scallop(config *_conf);
	~scallop();

public:
	int process(const char *bam_file);
	int add_bundle(bundle *bd);

};

#endif
