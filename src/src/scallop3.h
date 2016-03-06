#ifndef __SCALLOP3_H__
#define __SCALLOP3_H__

#include "assembler.h"
#include "algebra.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;

// algorithm: identify subsetsum signal
class scallop3 : public assembler
{
public:
	scallop3(splice_graph &gr);
	virtual ~scallop3();

public:
	MEI e2i;		// edge map
	VE i2e;			// edge map
	MEV mev;		// super edges
	VVD ns;			// null space from nodes

public:
	int assemble();
	int print();

private:
	int init_super_edges();
	int reconstruct_splice_graph();
	bool decompose_trivial_vertex(int x);
	int build_null_space();
};

#endif
