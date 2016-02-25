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
	int process(const string &file);

	int assemble_bam(const string &file);
	int assemble_gtf(const string &file);
	int assemble_example(const string &file);
};

#endif
