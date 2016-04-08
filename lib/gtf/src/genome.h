#ifndef __GTF_FILE_H__
#define __GTF_FILE_H__

#include <string>
#include <map>
#include "gene.h"

using namespace std;

class genome
{
public:
	genome(const string &file);
	~genome();

public:
	vector<gene> genes;
	map<string, int> s2i;

public:
	int read(const string &file);
	int write(const string &file) const;
	int build_index();
	const gene* get_gene(string name) const;
};

#endif
