#ifndef __SIMULATOR_H__
#define __SIMULATOR_H__

#include "gene.h"
#include "item.h"

using namespace std;

class simulator
{
public:
	simulator(int _num_exons, int _num_transcripts, int _max_length, int _max_expression);

public:
	int num_exons;
	int num_transcripts;
	int max_expression;
	int max_length;
	
public:
	int simulate_transcript(const string &tid, const string &gid, transcript &t);
	int simulate_gene(const string &gid, gene &g);
};

#endif
