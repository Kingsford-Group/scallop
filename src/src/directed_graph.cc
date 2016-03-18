#include "directed_graph.h"
#include "draw.h"

#include <climits>
#include <cstdio>
#include <cassert>
#include <tuple>
#include <algorithm>

using namespace std;

directed_graph::directed_graph()
{}

directed_graph::directed_graph(const graph_base &gr)
{
	clear();
	for(int i = 0; i < gr.num_vertices(); i++) add_vertex();

	PEE p = gr.edges();
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		add_edge((*it)->source(), (*it)->target());
	}
}

directed_graph::~directed_graph()
{
}


PEB directed_graph::edge(int s, int t) const
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

int directed_graph::exchange(int x, int y, int z)
{
	assert(x >= 0 && x < num_vertices());
	assert(y >= 0 && y < num_vertices());
	assert(z >= 0 && z < num_vertices());

	edge_iterator it1, it2;
	set<edge_descriptor> m;
	for(tie(it1, it2) = out_edges(x); it1 != it2; it1++)
	{
		if(check_path((*it1)->target(), y) == false) continue;
		m.insert(*it1);
	}
	for(tie(it1, it2) = in_edges(y); it1 != it2; it1++)
	{
		if(check_path(x, (*it1)->source()) == false) continue;
		m.insert(*it1);
	}
	for(tie(it1, it2) = out_edges(y); it1 != it2; it1++)
	{
		if(check_path((*it1)->target(), z) == false) continue;
		m.insert(*it1);
	}
	for(tie(it1, it2) = in_edges(z); it1 != it2; it1++)
	{
		if(check_path(y, (*it1)->source()) == false) continue;
		m.insert(*it1);
	}

	for(set<edge_descriptor>::iterator it = m.begin(); it != m.end(); it++)
	{
		int s = (*it)->source();
		int t = (*it)->target();
		if(s == x && t == y) move_edge(*it, y, z);
		else if(s == x) move_edge(*it, y, t);
		else if(t == y) move_edge(*it, s, z);
		else if(s == y && t == z) move_edge(*it, x, y);
		else if(s == y) move_edge(*it, x, t);
		else if(t == z) move_edge(*it, s, y);
		else assert(false);
	}
	return 0;
}

int directed_graph::remove_edge(edge_descriptor e)
{
	return graph_base::remove_edge(e);
}

int directed_graph::remove_edge(int s, int t)
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
		graph_base::remove_edge(v[i]);
	}
	return 0;
}

int directed_graph::in_degree(int v) const
{
	return vv[v]->in_degree();
}

int directed_graph::out_degree(int v) const
{
	return vv[v]->out_degree();
}

PEE directed_graph::edges() const
{
	return graph_base::edges();
}

PEE directed_graph::in_edges(int v) const
{
	return vv[v]->in_edges();
}

PEE directed_graph::out_edges(int v) const
{
	return vv[v]->out_edges();
}

vector<edge_descriptor> directed_graph::edges(int s, int t) const
{
	vector<edge_descriptor> v;
	PEE p = out_edges(s);
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		if((*it)->target() != t) continue;
		v.push_back(*it);
	}
	return v;
}

set<int> directed_graph::adjacent_vertices(int v) const
{
	return vv[v]->adjacent_vertices();
}

int directed_graph::move_edge(edge_base *e, int x, int y)
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

int directed_graph::bfs_reverse(int t, vector<int> &v, vector<int> &b) const
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
			int y = (*it1)->source();
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

int directed_graph::bfs_reverse(int t, vector<int> &v) const
{
	vector<int> b;
	return bfs_reverse(t, v, b);
}

bool directed_graph::check_path(int s, int t) const
{
	return graph_base::check_path(s, t);
}

bool directed_graph::check_path(edge_descriptor ex, edge_descriptor ey) const
{
	int s = ex->target();
	int t = ey->source();
	return check_path(s, t);
}

bool directed_graph::compute_shortest_path(int x, int y, vector<int> &p) const
{
	return graph_base::compute_shortest_path(x, y, p);
}

bool directed_graph::compute_shortest_path(edge_descriptor ex, edge_descriptor ey, vector<int> &p) const
{
	int s = ex->target();
	int t = ey->source();
	return graph_base::compute_shortest_path(s, t, p);
}

bool directed_graph::intersect(edge_descriptor ex, edge_descriptor ey) const
{
	int xs = ex->source();
	int xt = ex->target();
	int ys = ey->source();
	int yt = ey->target();
	if(xs == ys || xs == yt) return false;
	if(xt == ys || xt == yt) return false;
	if(check_path(xs, ys) == true)
	{
		if(check_path(ys, xt) == false) return false;
		if(check_path(xt, yt) == false) return false;
		return true;
	}
	if(check_path(ys, xs) == true)
	{
		if(check_path(xs, yt) == false) return false;
		if(check_path(yt, xt) == false) return false;
		return true;
	}
	return false;
}

vector<int> directed_graph::topological_sort() const
{
	vector<int> v;

	vector<int> vd;
	for(int i = 0; i < num_vertices(); i++)
	{
		int d = in_degree(i);
		vd.push_back(d);
	}

	vector<bool> vb;
	vb.assign(num_vertices(), false);
	for(int i = 0; i < num_vertices(); i++)
	{
		int k = -1;
		int mind = INT_MAX;
		for(int j = 0; j < num_vertices(); j++)
		{
			if(vb[j] == true) continue;
			if(vd[j] < mind)
			{
				mind = vd[j];
				k = j;
			}
		}
		assert(k != -1);
		v.push_back(k);
		vb[k] = true;

		edge_iterator it1, it2;
		for(tie(it1, it2) = out_edges(k); it1 != it2; it1++)
		{
			int t = (*it1)->target();
			vd[t]--;
		}
	}
	return v;
}

int directed_graph::compute_in_partner(int x) const
{
	vector<int> v = topological_sort();
	assert(v.size() == num_vertices());
	vector<int> order;
	order.assign(num_vertices(), -1);
	for(int i = 0; i < v.size(); i++)
	{
		order[v[i]] = i;
	}

	set<edge_descriptor> se;
	set<int> sv;
	edge_iterator it1, it2;
	for(tie(it1, it2) = in_edges(x); it1 != it2; it1++)
	{
		assert((*it1)->target() == x);
		int s = (*it1)->source();
		if(sv.find(s) == sv.end()) sv.insert(s);
		if(se.find(*it1) == se.end()) se.insert(*it1);
	}

	while(sv.size() >= 2)
	{
		int k = -1;
		for(set<int>::iterator it = sv.begin(); it != sv.end(); it++)
		{
			if(k == -1 || order[*it] > order[k])
			{
				k = *it;
			}
		}
		assert(k != -1);

		for(tie(it1, it2) = out_edges(k); it1 != it2; it1++)
		{
			if(se.find(*it1) == se.end()) return -1;
		}

		for(tie(it1, it2) = in_edges(k); it1 != it2; it1++)
		{
			assert(se.find(*it1) == se.end());
			se.insert(*it1);
			int s = (*it1)->source();
			if(sv.find(s) == sv.end()) sv.insert(s);
		}
		sv.erase(k);
	}
	if(sv.size() == 0) return -1;
	else return *(sv.begin());
}

int directed_graph::compute_out_partner(int x) const
{
	vector<int> v = topological_sort();
	vector<int> order;
	order.assign(num_vertices(), -1);
	for(int i = 0; i < v.size(); i++)
	{
		order[v[i]] = i;
	}

	set<edge_descriptor> se;
	set<int> sv;
	edge_iterator it1, it2;
	for(tie(it1, it2) = out_edges(x); it1 != it2; it1++)
	{
		assert((*it1)->source() == x);
		int t = (*it1)->target();
		if(sv.find(t) == sv.end()) sv.insert(t);
		if(se.find(*it1) == se.end()) se.insert(*it1);
	}

	while(sv.size() >= 2)
	{
		int k = -1;
		for(set<int>::iterator it = sv.begin(); it != sv.end(); it++)
		{
			if(k == -1 || order[*it] < order[k])
			{
				k = *it;
			}
		}
		assert(k != -1);

		for(tie(it1, it2) = in_edges(k); it1 != it2; it1++)
		{
			if(se.find(*it1) == se.end()) return -1;
		}

		for(tie(it1, it2) = out_edges(k); it1 != it2; it1++)
		{
			assert(se.find(*it1) == se.end());
			se.insert(*it1);
			int t = (*it1)->target();
			if(sv.find(t) == sv.end()) sv.insert(t);
		}
		sv.erase(k);
	}

	if(sv.size() == 0) return -1;
	else return *(sv.begin());
}

int directed_graph::draw(const string &file, const MIS &mis, const MES &mes, double len) const
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

	vector<int> tp = topological_sort();
	for(int ii = 0; ii < tp.size(); ii++)
	{
		int i = tp[ii];
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
			//assert(i < j);

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

			string line = "line width = 0.02cm,";
			if(cnt >= 2) line = "thin, double,";
			fout<<"\\draw[->,"<< line.c_str() <<"\\colx, bend right = "<< bend <<"] ("<<sx<<") to node[draw, fill=white] {";
			fout<< s.c_str() <<"} ("<<sy<<");\n";
		}
	}

	draw_footer(fout);

	fout.close();
	return 0;
}

int directed_graph::print() const
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
