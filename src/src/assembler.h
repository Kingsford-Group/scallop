#ifndef __ASSEMBLER_H__
#define __ASSEMBLER_H__

#include <fstream>
#include <string>

using namespace std;

class assembler
{
public:
	assembler();
	~assembler();

public:
	int process();

	int assemble_sgr(const string &file);
	int assemble_bam(const string &file);
	int assemble_gtf(const string &file);
};

#endif
