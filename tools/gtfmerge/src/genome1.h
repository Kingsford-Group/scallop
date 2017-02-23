#ifndef __GENOME1_H__
#define __GENOME1_H__

#include "genome.h"

using namespace std;

typedef pair<PI32, set<int> > PPIS;
typedef map<PI32, set<int> > MPIS;

class genome1: public genome
{
public:
	genome1();
	genome1(const string &file);

public:
	vector<transcript> transcripts;
	MPIS intron_index;	

public:
	int build();
	int print(int index);
	int query(const transcript &t, const set<int> &fb);

private:
	int collect_multiexon_transcripts();
	int build_intron_index();
};

#endif
