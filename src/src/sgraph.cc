#include "sgraph.h"

sgraph::sgraph()
{}

sgraph::~sgraph()
{}

int sgraph::add_node(double ave_abd, double dev_abd)
{
	add_vertex(gr);
	ave.push_back(ave_abd);
	dev.push_back(dev_abd);
	return 0;
}

int sgraph::add_arc(int s, int t, double wt)
{
	assert(s >= 0 && s < num_vertices(gr));
	assert(t >= 0 && t < num_vertices(gr));
	pair<edge_descriptor, bool> p = add_edge(s, t, gr);
	assert(p.second == false);
	MED::iterator it = ewt.find(p.first);
	assert(it == ewt.end());
	it->second = wt;
	return 0;
}
