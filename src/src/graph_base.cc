#include "graph_base.h"
#include "draw.h"

#include <climits>
#include <cstdio>
#include <cassert>
#include <tuple>
#include <algorithm>

using namespace std;

graph_base::graph_base()
{}

graph_base::~graph_base()
{
	for(int i = 0; i < vv.size(); i++) delete vv[i];
	for(edge_iterator it = se.begin(); it != se.end(); it++) delete (*it);
}

int graph_base::copy(const graph_base &gr)
{
	clear();
	for(int i = 0; i < gr.num_vertices(); i++) add_vertex();

	PEE p = gr.edges();
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		add_edge((*it)->source(), (*it)->target());
	}
	return 0;
}

int graph_base::add_vertex()
{
	vertex_base *v = new vertex_base();
	vv.push_back(v);
	return 0;
}

int graph_base::clear_vertex(int x)
{
	vector<edge_base*> v;
	PEE pi = vv[x]->in_edges();
	PEE po = vv[x]->out_edges();
	for(edge_iterator it = pi.first; it != pi.second; it++)
	{
		v.push_back(*it);
	}
	for(edge_iterator it = po.first; it != po.second; it++)
	{
		v.push_back(*it);
	}
	for(int i = 0; i < v.size(); i++)
	{
		remove_edge(v[i]);
	}
	return 0;
}

int graph_base::clear()
{
	for(int i = 0; i < vv.size(); i++) delete vv[i];
	for(edge_iterator it = se.begin(); it != se.end(); it++)
	{
		delete (*it);
	}
	vv.clear();
	se.clear();
	return 0;
}

int graph_base::degree(int v) const
{
	return vv[v]->degree();
}

PEB graph_base::edge(int s, int t) 
{
	assert(s >= 0 && s < vv.size());
	assert(t >= 0 && t < vv.size());
	PEE p = vv[s]->out_edges();
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		assert((*it)->source() == s);
		int x = (*it)->target();
		if(x != t) continue;
		return PEB(*it, true);
	}
	return PEB(null_edge, false);
}

vector<edge_descriptor> graph_base::edges(int s, int t)
{
	vector<edge_descriptor> v;
	PEE p = out_edges(s);
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		assert((*it)->source() == s);
		if((*it)->target() != t) continue;
		v.push_back(*it);
	}
	return v;
}

PEE graph_base::edges() const
{
	return PEE(se.begin(), se.end());
}

set<int> graph_base::adjacent_vertices(int s)
{
	set<int> ss;
	PEE p = out_edges(s);
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		assert((*it)->source() == s);
		int t = (*it)->target();
		if(ss.find(t) == ss.end()) ss.insert(t);
	}
	return ss;
}

size_t graph_base::num_vertices() const
{
	return vv.size();
}

size_t graph_base::num_edges() const
{
	return se.size();
}

int graph_base::bfs(int s, vector<int> &v, vector<int> &b)
{
	v.assign(num_vertices(), -1);
	b.assign(num_vertices(), -1);
	vector<bool> closed;
	closed.resize(num_vertices(), false);
	vector<int> open;
	open.push_back(s);
	v[s] = 0;
	int p = 0;

	while(p < open.size())
	{
		int x = open[p];
		assert(v[x] >= 0);

		p++;
		if(closed[x] == true) continue;
		closed[x] = true;

		edge_iterator it1, it2;
		for(tie(it1, it2) = out_edges(x); it1 != it2; it1++)
		{
			int y = (*it1)->target();
			if(v[y] == -1) 
			{
				v[y] = 1 + v[x];
				b[y] = x;
			}
			assert(v[y] <= 1 + v[x]);
			if(closed[y] == true) continue;
			open.push_back(y);
		}
	}
	return 0;
}

int graph_base::bfs(int t, vector<int> &v) 
{
	vector<int> b;
	return bfs(t, v, b);
}

bool graph_base::check_path(int s, int t) 
{
	vector<int> v;
	bfs(s, v);
	if(v[t] == -1) return false;
	else return true;
}

bool graph_base::compute_shortest_path(int s, int t, vector<int> &p) 
{
	p.clear();
	vector<int> v;
	vector<int> b;
	bfs(s, v, b);
	if(v[t] == -1) return false;
	int x = t;
	while(true)
	{
		p.push_back(x);
		if(x == s) break;
		x = b[x];
	}
	reverse(p.begin(), p.end());
	return true;
}

bool graph_base::check_nested() 
{
	PEE p = edges();
	for(edge_iterator i = p.first; i != p.second; i++)
	{
		edge_iterator j = i;
		j++;
		for(; j != p.second; j++)
		{
			bool b = intersect(*i, *j);
			if(b == true) return false;
		}
	}
	return true;
}

int graph_base::print() const
{
	printf("total %lu vertices, %lu edges\n", vv.size(), se.size());
	for(int i = 0; i < vv.size(); i++)
	{
		printf("vertex %d: ", i);
		vv[i]->print();
	}

	for(edge_iterator it = se.begin(); it != se.end(); it++)
	{
		(*it)->print();
	}
	return 0;
}

