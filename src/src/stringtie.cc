#include "stringtie.h"

stringtie::stringtie(string &s, splice_graph &gr)
	: assembler(s, gr)
{}

stringtie::~stringtie()
{}

int stringtie::assemble()
{
	smooth_weights();
	greedy();
	printf("%s solution %lu paths\n", name.c_str(), paths.size());
	return 0;
}

int stringtie::greedy()
{
	while(true)
	{
		path p = compute_maximum_forward_path();
		decrease_path(p);
		if(p.v.size() <= 1) break;
		paths.push_back(p);
		if(p.abd < 1) break;
	}
	return 0;
}

