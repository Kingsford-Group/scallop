#ifndef __GTF_FILE_H__
#define __GTF_FILE_H__

#include <string>
#include <map>
#include "gene.h"

using namespace std;

class genome
{
public:
	genome();
	genome(const string &file);
	~genome();

public:
	vector<gene> genes;
	map<string, int> s2i;

public:
	// read and write
	int read(const string &file);
	int write(const string &file) const;

	// modify
	int add_gene(const gene &g);
	int sort();
	int build_index();
	int assign_RPKM(int reads);

	// fetch information
	const gene* get_gene(string name) const;
	const gene* locate_gene(const string &chr, const PI32 &p) const;

};

#endif
