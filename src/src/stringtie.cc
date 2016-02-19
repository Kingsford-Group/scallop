#include "stringtie.h"

stringtie::stringtie(splice_graph &gr)
	: assembler(gr)
{}

stringtie::~stringtie()
{}

int stringtie::assemble()
{
	update_weights();
	greedy();
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

