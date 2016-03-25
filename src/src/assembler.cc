#include "assembler.h"
#include "smoother.h"


assembler::assembler(const string &s, const splice_graph &g)
	: gr(g), name(s)
{}

assembler::~assembler()
{}

int assembler::smooth_weights()
{
	smoother lp(gr);
	lp.solve();
	return 0;
}

int assembler::print(const string &prefix) const
{
	// print paths
	double w = 0.0;
	for(int i = 0; i < paths.size(); i++)
	{
		printf("%s ", prefix.c_str());
		paths[i].print(i);
		w += paths[i].abd;
	}
	printf("%s summary: %ld paths, %.2lf total abundance\n", prefix.c_str(), paths.size(), w);

	return 0;
}
