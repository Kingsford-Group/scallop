#include "graph_base.h"
#include "draw.h"

#include <cstdio>
#include <cassert>
#include <tuple>
#include <algorithm>

using namespace std;

graph_base::graph_base()
{}

graph_base::graph_base(const graph_base &gr)
{
	clear();
	for(int i = 0; i < gr.num_vertices(); i++) add_vertex();

	PEE p = gr.edges();
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		add_edge((*it)->source(), (*it)->target());
	}
}

graph_base::~graph_base()
{
	for(int i = 0; i < vv.size(); i++) delete vv[i];
	for(edge_iterator it = se.begin(); it != se.end(); it++) delete (*it);
}

int graph_base::add_vertex()
{
	vertex_base *v = new vertex_base();
	vv.push_back(v);
	return 0;
}

edge_descriptor graph_base::add_edge(int s, int t)
{
	assert(s >= 0 && s < vv.size());
	assert(t >= 0 && t < vv.size());
	edge_base *e = new edge_base(s, t);
	assert(se.find(e) == se.end());
	se.insert(e);
	vv[s]->add_out_edge(e);
	vv[t]->add_in_edge(e);
	return e;
}

PEB graph_base::edge(int s, int t) const
{
	assert(s >= 0 && s < vv.size());
	assert(t >= 0 && t < vv.size());
	PEE p = vv[s]->out_edges();
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		int x = (*it)->target();
		if(x != t) continue;
		return PEB(*it, true);
	}
	return PEB(null_edge, false);
}

int graph_base::remove_edge(edge_base *e)
{
	if(se.find(e) == se.end()) return -1;
	vv[e->source()]->remove_out_edge(e);
	vv[e->target()]->remove_in_edge(e);
	delete e;
	se.erase(e);
	return 0;
}

int graph_base::remove_edge(int s, int t)
{
	vector<edge_base*> v;
	PEE p = vv[s]->out_edges();
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		if((*it)->target() != t) continue;
		v.push_back(*it);
	}
	for(int i = 0; i < v.size(); i++)
	{
		remove_edge(v[i]);
	}
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

int graph_base::in_degree(int v) const
{
	return vv[v]->in_degree();
}

int graph_base::out_degree(int v) const
{
	return vv[v]->out_degree();
}

PEE graph_base::in_edges(int v) const
{
	return vv[v]->in_edges();
}

PEE graph_base::out_edges(int v) const
{
	return vv[v]->out_edges();
}

PEE graph_base::edges() const
{
	return PEE(se.begin(), se.end());
}


size_t graph_base::num_vertices() const
{
	return vv.size();
}

size_t graph_base::num_edges() const
{
	return se.size();
}

set<int> graph_base::adjacent_vertices(int v) const
{
	return vv[v]->adjacent_vertices();
}

int graph_base::move_edge(edge_base *e, int x, int y)
{
	int s = e->source();
	int t = e->target();

	if(s != x)
	{
		vv[s]->remove_out_edge(e);
		vv[x]->add_out_edge(e);
	}

	if(t != y)
	{
		vv[t]->remove_in_edge(e);
		vv[y]->add_in_edge(e);
	}

	e->move(x, y);
	return 0;
}

int graph_base::bfs(int s, vector<int> &v, vector<int> &b) const
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


int graph_base::bfs_reverse(int t, vector<int> &v, vector<int> &b) const
{
	v.assign(num_vertices(), -1);
	b.assign(num_vertices(), -1);
	vector<bool> closed;
	closed.resize(num_vertices(), false);
	vector<int> open;
	open.push_back(t);
	v[t] = 0;
	int p = 0;

	while(p < open.size())
	{
		int x = open[p];
		assert(v[x] >= 0);

		p++;
		if(closed[x] == true) continue;
		closed[x] = true;

		edge_iterator it1, it2;
		for(tie(it1, it2) = in_edges(x); it1 != it2; it1++)
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

int graph_base::bfs(int t, vector<int> &v) const
{
	vector<int> b;
	return bfs(t, v, b);
}

int graph_base::bfs_reverse(int t, vector<int> &v) const
{
	vector<int> b;
	return bfs_reverse(t, v, b);
}

bool graph_base::check_directed_path(int s, int t) const
{
	vector<int> v;
	bfs(s, v);
	if(v[t] == -1) return false;
	else return true;
}

bool graph_base::check_directed_path(edge_descriptor ex, edge_descriptor ey) const
{
	int s = ex->target();
	int t = ey->source();
	return check_directed_path(s, t);
}

bool graph_base::compute_shortest_path(int s, int t, vector<int> &p) const
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

bool graph_base::compute_shortest_path(edge_descriptor ex, edge_descriptor ey, vector<int> &p) const
{
	int s = ex->target();
	int t = ey->source();
	return compute_shortest_path(s, t, p);
}

bool graph_base::intersect(edge_descriptor ex, edge_descriptor ey) const
{
	int xs = ex->source();
	int xt = ex->target();
	int ys = ey->source();
	int yt = ey->target();
	if(xs == ys || xs == yt) return false;
	if(xt == ys || xt == yt) return false;
	if(check_directed_path(xs, ys) == true)
	{
		if(check_directed_path(ys, xt) == false) return false;
		if(check_directed_path(xt, yt) == false) return false;
		return true;
	}
	if(check_directed_path(ys, xs) == true)
	{
		if(check_directed_path(xs, yt) == false) return false;
		if(check_directed_path(yt, xt) == false) return false;
		return true;
	}
	return false;
}

bool graph_base::check_nested() const
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

int graph_base::draw(const string &file, const MIS &mis, const MES &mes, double len) const
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
		string s = "";
		MIS::const_iterator it = mis.find(i);
		if(it != mis.end()) s = it->second;
		fout.precision(0);
		fout<<fixed;
		fout<<"\\node[mycircle, \\colx, draw, label = below:{";
		fout<<s.c_str() <<"}] ("<<sx<<") at ("<<pos<<" *\\len, 0.0) {"<<i<<"};\n";
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

			string s;
			char buf[1024];
			edge_iterator oi1, oi2;
			int cnt = 0;
			for(tie(oi1, oi2) = out_edges(i); oi1 != oi2; oi1++)
			{
				if((*oi1)->target() != j) continue;
				cnt++;
				MES::const_iterator it = mes.find(*oi1);
				if(it == mes.end()) continue;
				if(cnt == 1) s.append(it->second);
				else s.append("," + it->second);
			}

			sprintf(sx, "s%d", i);
			sprintf(sy, "s%d", j);

			double bend = -40;
			if(i + 1 == j) bend = 0;

			string line = "";
			if(cnt == 1) line = "line width = 0.02cm,";
			else line = "thin, double,";
			fout<<"\\draw[->,"<< line.c_str() <<"\\colx, bend right = "<< bend <<"] ("<<sx<<") to node {";
			fout<< s.c_str() <<"} ("<<sy<<");\n";
		}
	}

	draw_footer(fout);

	fout.close();
	return 0;
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

int graph_base::test()
{
	graph_base gr;
	gr.add_vertex();
	gr.add_vertex();
	gr.add_vertex();
	gr.add_vertex();
	gr.add_vertex();

	gr.add_edge(0, 1);
	gr.add_edge(0, 1);
	gr.add_edge(0, 1);
	gr.add_edge(1, 2);
	gr.add_edge(2, 3);
	gr.add_edge(2, 4);
	gr.add_edge(3, 4);

	gr.print();

	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		gr.remove_edge(*it1);
	}

	gr.print();

	return 0;
}
