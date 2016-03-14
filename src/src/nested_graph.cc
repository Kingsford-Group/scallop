#include "nested_graph.h"

nested_graph::nested_graph()
{}

nested_graph::nested_graph(const graph_base &gr)
{
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		add_vertex();
	}

	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		add_edge(s, t);
	}
}

nested_graph::~nested_graph()
{}

edge_descriptor nested_graph::add_edge(int s, int t)
{
	edge_descriptor e = graph_base::add_edge(s, t);
	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(); it1 != it2; it1++)
	{
		edge_descriptor ee = *it1;
		if(e == ee) continue;
		bool b = intersect(e, ee);
		if(b == false) continue;
		int ss = e->source() < ee->source() ? e->source() : ee->source();
		int tt = e->target() > ee->target() ? e->target() : ee->target();
		remove_edge(e);
		remove_edge(ee);
		add_edge(ss, tt);
		return null_edge;
	}
	return e;
}

bool nested_graph::intersect(edge_descriptor &ex, edge_descriptor &ey) const
{
	int xs = ex->source();
	int xt = ex->target();
	int ys = ey->source();
	int yt = ey->target();
	assert(xs < xt);
	assert(ys < yt);

	if(xs == ys) return false;
	if(xs > ys) return intersect(ey, ex);
	if(ys >= xt) return false;
	if(yt <= xt) return false;
	// TODO, right now this is a over-strong condition
	return true;
}

