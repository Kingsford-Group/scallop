#ifndef __METAASSEMBLER_H__
#define __METAASSEMBLER_H__

#include "assembler.h"

class metaassembler
{
public:
	vector<splice_graph> grlist1;
	vector<splice_graph> grlist2;
	vector<splice_graph> grlist3;
	vector<hyper_set> hslist1;
	vector<hyper_set> hslist2;
	vector<hyper_set> hslist3;

public:
	int library_type1;
	int library_type2;

public:
	int preassemble();
	int assemble();
	int postassemble();

	int print();
};

#endif
