#ifndef __ASSEMBLER_H__
#define __ASSEMBLER_H__

#include <fstream>
#include <string>
#include "bundle_base.h"
#include "bundle.h"

using namespace std;

class assembler
{
public:
	assembler();
	~assembler();

private:
	samFile *sfn;
	bam_hdr_t *hdr;
	bam1_t *b1t;
	
	ofstream fout;
	int index;

	vector<bundle_base> vbb;

	bool terminate;

public:
	int assemble();

private:
	int add_hit(const hit &ht);
	int truncate(const hit &ht);
	int process(const bundle_base &bb);
	int compare(splice_graph &gr, const string &ref, const string &tex);
};

#endif
