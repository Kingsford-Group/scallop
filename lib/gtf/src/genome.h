/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

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
	virtual ~genome();

public:
	vector<gene> genes;
	map<string, int> g2i;

public:
	// read and write
	int read(const string &file);
	int write(const string &file) const;

	// modify
	int add_gene(const gene &g);
	int sort();
	int build_index();
	int assign_RPKM(double factor);
	int assign_TPM_by_RPKM();
	int assign_TPM_by_FPKM();

	// filter
	int filter_single_exon_transcripts();
	int filter_low_coverage_transcripts(double min_coverage);

	// fetch information
	const gene* get_gene(string name) const;
	const gene* locate_gene(const string &chr, const PI32 &p) const;
	vector<transcript> collect_transcripts() const;
};

#endif
