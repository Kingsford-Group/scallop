#include "assembler.h"
#include "smoother.h"
#include "fibonacci_heap.h"

#include <iomanip>
#include <cfloat>

assembler::assembler(const splice_graph &g)
	: gr(g)
{}

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

path assembler::compute_maximum_forward_path()
{
	path p;
	vector<double> table;		// dynamic programming table
	vector<int> back;			// back pointers
	table.resize(gr.num_vertices(), 0);
	back.resize(gr.num_vertices(), -1);
	table[0] = DBL_MAX;
	int n = gr.num_vertices();
	for(int i = 1; i < n; i++)
	{
		double abd = gr.get_vertex_weight(i);
		if(i == n - 1) abd = DBL_MAX;

		double max_abd = 0;
		int max_idx = -1;
		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.in_edges(i); it1 != it2; it1++)
		{
			int s = (*it1)->source();
			int t = (*it1)->target();
			assert(t == i);
			assert(s < i);
			//if(s >= i) continue;
			double xw = gr.get_edge_weight(*it1);
			double ww = xw < table[s] ? xw : table[s];
			if(ww >= max_abd)
			{
				max_abd = ww;
				max_idx = s;
			}
		}

		back[i] = max_idx;
		table[i] = max_abd < abd ? max_abd : abd;
	}

	if(table[n - 1] <= 0) return p;

	p.abd = table[n - 1];
	p.v.clear();

	int b = n - 1;
	while(true)
	{
		p.v.push_back(b);
		if(b == 0) break;
		b = back[b];
		assert(b != -1);
	}
	reverse(p.v.begin(), p.v.end());

	return p;
}

path assembler::compute_maximum_path()
{
	path p;
	fibonacci_heap f;						// fibonacci heap
	vector<handle_t> handles;		// handles for nodes
	vector<bool> reached;			// whether each node is reached
	vector<bool> visited;			// whether each node is visited
	vector<int> backptr;			// backward pointers

	handles.resize(gr.num_vertices());
	reached.resize(gr.num_vertices(), false);
	visited.resize(gr.num_vertices(), false);
	backptr.resize(gr.num_vertices(), -1);

	handle_t h = f.push(fnode(0, DBL_MAX));
	handles[0] = h;
	reached[0] = true;
	backptr[0] = -1;

	while(f.empty() == false)
	{
		fnode fx = f.top();
		visited[fx.v] = true;

		if(fx.v == gr.num_vertices() - 1) break;
		
		f.pop();
		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.out_edges(fx.v); it1 != it2; it1++)
		{
			int y = (*it1)->target();
			if(visited[y] == true) continue;

			double xw = gr.get_edge_weight(*it1);
			double ww = xw < fx.w ? xw : fx.w;
			if(reached[y] == false) 
			{
				h = f.push(fnode(y, ww));
				handles[y] = h;
				reached[y] = true;
				backptr[y] = fx.v;
			}
			else
			{
				double yw = (*handles[y]).w;
				if(ww > yw) 
				{
					f.increase(handles[y], fnode(y, ww));
					backptr[y] = fx.v;
				}
			}
		}
	}

	if(f.empty() == true) return p;

	p.abd = f.top().w;
	p.v.clear();
	int b = gr.num_vertices() - 1;
	while(true)
	{
		p.v.push_back(b);
		if(b == 0) break;
		b = backptr[b];
		assert(b != -1);
	}
	reverse(p.v.begin(), p.v.end());

	return p;
}

int assembler::decrease_path(const path &p)
{
	if(p.v.size() < 2) return 0;
	for(int i = 0; i < p.v.size() - 1; i++)
	{
		PEB e = gr.edge(p.v[i], p.v[i + 1]);
		assert(e.second == true);
		double w0 = gr.get_edge_weight(e.first);
		double w1 = w0 - p.abd;
		assert(w1 >= -0.000001);
		if(w1 <= 0) w1 = 0;
		gr.set_edge_weight(e.first, w1);
	}
	return 0;
}

int assembler::increase_path(const path &p)
{
	if(p.v.size() < 2) return 0;
	for(int i = 0; i < p.v.size() - 1; i++)
	{
		PEB e = gr.edge(p.v[i], p.v[i + 1]);
		assert(e.second == true);
		double w0 = gr.get_edge_weight(e.first);
		double w1 = w0 + p.abd;
		gr.set_edge_weight(e.first, w1);
	}
	return 0;
}

int assembler::add_backward_path(const path &p)
{
	if(p.v.size() < 2) return 0;
	for(int i = 0; i < p.v.size() - 1; i++)
	{
		edge_descriptor e = gr.add_edge(p.v[i + 1], p.v[i]);
		gr.set_edge_weight(e, p.abd);
	}
	return 0;
}

int assembler::remove_backward_path(const path &p)
{
	if(p.v.size() < 2) return 0;
	for(int i = 0; i < p.v.size() - 1; i++)
	{
		gr.edge(p.v[i + 1], p.v[i]);
	}
	return 0;
}

double assembler::compute_bottleneck_weight(const path &p)
{
	double ww = DBL_MAX;
	for(int i = 0; i < p.v.size() - 1; i++)
	{
		PEB e = gr.edge(p.v[i], p.v[i + 1]);
		assert(e.second == true);
		double w = gr.get_edge_weight(e.first);
		if(w < ww) ww = w;
	}
	return ww;
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
