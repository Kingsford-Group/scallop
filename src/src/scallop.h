#ifndef __SCALLOP_H__
#define __SCALLOP_H__

#include "bundle.h"
#include <fstream>

using namespace std;

class scallop
{
public:
	scallop();
	~scallop();

public:
	int process(const string &file);
	int solve(const char *bam_file);
	int solve_bundle(const bundle &bd, int index, ofstream &fout1, ofstream &fout2);
};

#endif
