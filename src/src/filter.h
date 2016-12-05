#ifndef __FILTER_H__
#define __FILTER_H__

#include "gene.h"

class filter
{
public:
	filter(const vector<transcript> &v);

public:
	vector<transcript> trs;

public:
	int join();

private:
	bool join_transcripts();
	int locate_next_transcript(int t);
};

bool transcript_cmp(const transcript &x, const transcript &y);

#endif
