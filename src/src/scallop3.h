#ifndef __SCALLOP3_H__
#define __SCALLOP3_H__

#include "assembler.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;

// dynamic programming for nested graph
class scallop3 : public assembler
{
public:
	scallop3(splice_graph &gr);
	virtual ~scallop3();

public:
	MEV mev;

public:
	int assemble();
	bool decide_nested();

private:
	int print();
	int init_super_edges();
	int reconstruct_splice_graph();
	bool decompose_trivial_vertex(int x);
};

#endif
