#ifndef __MANAGER_H__
#define __MANAGER_H__

#include <fstream>

using namespace std;

class manager
{
public:
	manager();
	~manager();

public:
	ofstream standard_fout;
	ofstream stringtie_fout;
	ofstream scallop_fout;

public:
	int process(const string &file);

	int assemble_bam(const string &file);
	int assemble_gtf(const string &file);
	int assemble_example(const string &file);
};

#endif
