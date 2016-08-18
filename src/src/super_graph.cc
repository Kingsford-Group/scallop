#include "super_graph.h"

super_graph::super_graph()
{}

super_graph::~super_graph()
{}

int super_graph::build()
{
	build_undirected_graph();
	split_splice_graph();
	return 0;
}

int super_graph::build_undirected_graph()
{
	ug.clear();
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		ug.add_vertex();
	}

	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int s = e->source();
		int t = e->target();
		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;

		ug.add_edge(s, t);
	}

	return 0;
}

int super_graph::split_splice_graph()
{
	// TODO
	vector< set<int> > vv = ug.compute_connected_components();
	subgraphs.resize(vv.size());
	return 0;
}
