#ifndef __BUNDLE_H__
#define __BUNDLE_H__

#include <stdint.h>
#include <cstring>
#include <string>
#include <vector>

#include "config.h"
#include "hit.h"
#include "sam.h"

using namespace std;

class bundle
{
public:
	bundle(config *_conf);
	~bundle();

public:
	config *conf;				// config
	int32_t tid;				// chromosome ID
	string chrm;				// chromosome name
	int32_t lpos;				// the leftmost position on reference
	int32_t rpos;				// the rightmost position on reference
	vector<hit> hits;			// hits

public:
	int solve();
	int infer_splice_positions();
	int add_hit(bam_hdr_t *h, bam1_t *b);
	int print();
	int clear();
};

#endif
