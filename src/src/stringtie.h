#ifndef __STRINGTIE_H__
#define __STRINGTIE_H__

#include "assembler.h"

class stringtie : public assembler
{
public:
	stringtie(string &s, splice_graph &g);
	virtual ~stringtie();

public:
	int assemble();

private:
	int greedy();

private:
	double compute_bottleneck_weight(const path &p);
	path compute_maximum_forward_path();
	path compute_maximum_path();

	int decrease_path(const path &p);
	int increase_path(const path &p);
	int add_backward_path(const path &p);
	int remove_backward_path(const path &p);
};

#endif
