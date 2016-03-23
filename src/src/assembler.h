#ifndef __ASSEMBLER_H__
#define __ASSEMBLER_H__

#include "splice_graph.h"
#include "path.h"

class assembler
{
public:
	assembler(const string &s, const splice_graph &g);
	virtual ~assembler();

public:
	string name;					// name for this gene
	splice_graph gr;				// splice graph
	vector<path> paths;				// transcripts

public:
	virtual int assemble() = 0;
	int print(const string &prefix) const;

protected:
	int smooth_weights();

	double compute_bottleneck_weight(const path &p);
	path compute_maximum_forward_path();
	path compute_maximum_path();

	int decrease_path(const path &p);
	int increase_path(const path &p);
	int add_backward_path(const path &p);
	int remove_backward_path(const path &p);
};

#endif
