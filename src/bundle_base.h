/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __BUNDLE_BASE_H__
#define __BUNDLE_BASE_H__

#include <stdint.h>
#include <cstring>
#include <string>
#include <vector>

#include "hit.h"
#include "interval_map.h"
#include "transcript.h"
#include "config.h"

using namespace std;

class bundle_base
{
public:
	bundle_base();
	virtual ~bundle_base();

public:
	int32_t tid;					// chromosome ID
	string chrm;					// chromosome name
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	char strand;					// strandness
	vector<hit> hits;				// hits
	split_interval_map mmap;		// matched interval map
	split_interval_map imap;		// indel interval map
	vector<transcript> trsts;

public:
	int add_hit(const hit &ht);
	bool overlap(const hit &ht) const;
	int clear();
	vector<int> hits_on_transcripts(vector<uint32_t> &read_mapping);
};

#endif
