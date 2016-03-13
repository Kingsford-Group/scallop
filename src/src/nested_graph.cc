#include "nested_graph.h"
#include "draw.h"

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
		return add_edge(ss, tt);
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

int nested_graph::draw(const string &file, double len) const
{
	ofstream fout(file.c_str());
	if(fout.fail())
	{
		printf("open file %s error.\n", file.c_str());
		return 0;
	}

	draw_header(fout);

	fout<<"\\def\\len{"<<len<<"cm}\n";

	// draw vertices
	char sx[1024];
	char sy[1024];
	double pos = 0;
	for(int i = 0; i < num_vertices(); i++)
	{
		int d = degree(i);
		if(d == 0) continue;

		pos++;

		sprintf(sx, "s%d", i);
		fout.precision(0);
		fout<<fixed;
		fout<<"\\node[mycircle, \\colx, draw] ("<<sx<<") at ("<<pos<<" *\\len, 0.0) {"<<i<<"};\n";
	}

	// draw edges
	for(int i = 0; i < num_vertices(); i++)
	{
		set<int> ss = adjacent_vertices(i);
		for(set<int>::iterator it = ss.begin(); it != ss.end(); it++)
		{
			// TODO
			int j = (*it);
			assert(i < j);

			int cnt = 0;
			edge_iterator oi1, oi2;
			for(tie(oi1, oi2) = out_edges(i); oi1 != oi2; oi1++)
			{
				if((*oi1)->target() != j) continue;
				cnt++;
			}

			char buf[1024];
			sprintf(buf, "%d", cnt);

			sprintf(sx, "s%d", i);
			sprintf(sy, "s%d", j);

			double bend = -40;
			if(i + 1 == j) bend = 0;

			string line = "line width = 0.02cm,";

			fout<<"\\draw[->,"<< line.c_str() <<"\\colx, bend right = "<< bend <<"] ("<<sx<<") to node {";
			fout<< buf <<"} ("<<sy<<");\n";
		}
	}

	draw_footer(fout);

	fout.close();
	return 0;
}

