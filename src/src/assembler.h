#ifndef __ASSEMBLER_H__
#define __ASSEMBLER_H__

#include <fstream>
#include <string>
#include "bundle_base.h"

using namespace std;

class assembler
{
public:
	assembler();
	~assembler();

public:
	int process();

	int assemble_sgr(const string &file);
	int assemble_gtf(const string &file);
	int assemble_bam(const string &file);
	int process_bundle(bundle_base &bb, int &index, ofstream &fout);
};

#endif
