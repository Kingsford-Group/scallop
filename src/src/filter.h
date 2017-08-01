/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

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
	int join_single_exon_transcripts();
	int filter_length_coverage();
	int remove_nested_transcripts();
	int merge_single_exon_transcripts();
	int merge_single_exon_transcripts(vector<transcript> &trs0);
	int print();

private:
	bool join_transcripts();
	int locate_next_transcript(int t);
};

bool transcript_cmp(const transcript &x, const transcript &y);

#endif
