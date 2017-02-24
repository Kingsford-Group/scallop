#ifndef __GENOME1_H__
#define __GENOME1_H__

#include "genome.h"

using namespace std;

typedef pair<PI32, set<int> > PPIS;
typedef map<PI32, set<int> > MPIS;
typedef pair<int, int> PII;
typedef map<int, int> MII;

class genome1
{
public:
	genome1();
	genome1(const string &file);

public:
	vector<transcript> transcripts;
	MPIS intron_index;	

public:
	int build(const string &file);
	int build(const vector<transcript> &v);
	int build_union(const genome1 &gm);
	int build_intersection(const genome1 &gm, genome1 &out);
	int clear();
	int print(int index);
	int write(const string &file);
	int add_suffix(const string &p);

private:
	int build_multiexon_transcripts(const string &file);
	int build_intron_index();
	int query(const transcript &t, const set<int> &fb);
	int compare(const genome1 &gy, MII &x2y, MII &y2x);
};

#endif
