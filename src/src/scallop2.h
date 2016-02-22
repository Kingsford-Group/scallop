#ifndef __SCALLOP2_H__
#define __SCALLOP2_H__

#include "gurobi_c++.h"
#include "assembler.h"

class scallop2 : public assembler
{
public:
	scallop2(splice_graph &gr);
	~scallop2();

public:
	MED ewrt;						// edge weights

	double cabd;
	vector<path> cpaths;			// current paths
	MEI ecnt;						// number of paths covering e

	// for optimizing the abundance 
	typedef map<edge_descriptor, GRBLinExpr> MEG;
	typedef pair<edge_descriptor, GRBLinExpr> PEG;
	GRBEnv *env;					// GRBEnv
	GRBModel *model;				// GRBModel
	vector<GRBVar> vars;			// vars for paths
	GRBLinExpr obj;					// objective function
	MEG eexprs;						// sum of paths covering e

public:
	int assemble();

private:
	int init();
	int iterate();

	int assign_weights();
	int update_counters(const path &p);
	int update_lpsolver(const path &p);
	double optimize();
	int assign_abundance();
	int collect();
};

#endif
