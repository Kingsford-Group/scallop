/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "directed_graph.h"
#include "draw.h"

#include <climits>
#include <cstdio>
#include <cassert>
#include <tuple>
#include <algorithm>
#include <iostream>

using namespace std;

directed_graph::directed_graph()
{}

directed_graph::directed_graph(const directed_graph &gr)
{
	copy(gr);
}

directed_graph& directed_graph::operator=(const directed_graph &gr)
{
	copy(gr);
	return (*this);
}

directed_graph::~directed_graph()
{
}

edge_descriptor directed_graph::add_edge(int s, int t)
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

int directed_graph::remove_edge(edge_descriptor e)
{
	if(se.find(e) == se.end()) return -1;
	vv[e->source()]->remove_out_edge(e);
	vv[e->target()]->remove_in_edge(e);
	delete e;
	se.erase(e);
	return 0;
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

int directed_graph::rotate(int x, int y)
{
	if(check_path(y, x) == true) return rotate(y, x);

	set<edge_descriptor> se;
	int f = check_nest(x, y, se);
	assert(f >= 0);

	for(set<edge_descriptor>::iterator it = se.begin(); it != se.end(); it++)
	{
		int s = (*it)->source();
		int t = (*it)->target();

		if(s == x && t == y) continue;
		else if(s == x) move_edge(*it, t, y);
		else if(t == y) move_edge(*it, x, s);
		else move_edge(*it, t, s);
	}
	return 0;
}

int directed_graph::remove_edge(int s, int t)
{
	vector<edge_base*> v;
	PEEI p = vv[s]->out_edges();
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

int directed_graph::in_degree(int v) const
{
	return vv[v]->in_degree();
}

int directed_graph::out_degree(int v) const
{
	return vv[v]->out_degree();
}

PEEI directed_graph::in_edges(int v)
{
	return vv[v]->in_edges();
}

PEEI directed_graph::out_edges(int v)
{
	return vv[v]->out_edges();
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

bool directed_graph::bfs_reverse(const vector<int> &t, int s, const set<edge_descriptor> &fb)
{
	vector<int> open = t;
	set<int> closed(t.begin(), t.end());
	int p = 0;

	while(p < open.size())
	{
		int x = open[p];
		if(x == s) return true;
		p++;
		edge_iterator it1, it2;
		for(tie(it1, it2) = in_edges(x); it1 != it2; it1++)
		{
			if(fb.find(*it1) != fb.end()) continue;
			int y = (*it1)->source();
			if(closed.find(y) != closed.end()) continue;
			closed.insert(y);
			open.push_back(y);
		}
	}
	return false;
}

int directed_graph::bfs_reverse(int t, vector<int> &v, vector<int> &b)
{
	v.clear();
	b.clear();
	set<int> closed;
	vector<int> open;
	open.push_back(t);
	closed.insert(t);
	v.push_back(t);
	b.push_back(-1);
	int p = 0;

	while(p < open.size())
	{
		int x = open[p];
		p++;
		edge_iterator it1, it2;
		for(tie(it1, it2) = in_edges(x); it1 != it2; it1++)
		{
			int y = (*it1)->source();

			if(closed.find(y) != closed.end()) continue;
			closed.insert(y);
			open.push_back(y);

			v.push_back(y);
			b.push_back(x);
		}
	}
	return 0;
}

int directed_graph::bfs_reverse(int t, set<edge_descriptor> &ss)
{
	ss.clear();
	set<int> closed;
	vector<int> open;
	open.push_back(t);
	closed.insert(t);
	int p = 0;
	while(p < open.size())
	{
		int x = open[p];
		p++;
		edge_iterator it1, it2;
		for(tie(it1, it2) = in_edges(x); it1 != it2; it1++)
		{
			ss.insert(*it1);
			int y = (*it1)->source();
			if(closed.find(y) != closed.end()) continue;
			closed.insert(y);
			open.push_back(y);
		}
	}
	return 0;
}

int directed_graph::bfs_reverse(int t, vector<int> &v)
{
	vector<int> b;
	return bfs_reverse(t, v, b);
}

bool directed_graph::check_path(int s, int t)
{
	return graph_base::check_path(s, t);
}

bool directed_graph::check_path(edge_descriptor ex, edge_descriptor ey)
{
	int s = ex->target();
	int t = ey->source();
	return check_path(s, t);
}

bool directed_graph::compute_shortest_path(int x, int y, vector<int> &p)
{
	return graph_base::compute_shortest_path(x, y, p);
}

bool directed_graph::compute_shortest_path(edge_descriptor ex, edge_descriptor ey, vector<int> &p)
{
	int s = ex->target();
	int t = ey->source();
	return graph_base::compute_shortest_path(s, t, p);
}

bool directed_graph::intersect(edge_descriptor ex, edge_descriptor ey)
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

vector<int> directed_graph::topological_sort_reverse()
{
	vector<int> v;
	vector<int> q;
	vector<int> vd;
	for(int i = num_vertices() - 1; i >= 0; i--)
	{
		int d = out_degree(i);
		vd.push_back(d);
		if(d == 0) q.push_back(i);
	}

	int k = 0;

	while(k < q.size())
	{
		int x = q[k++];
		v.push_back(x);

		edge_iterator it1, it2;
		for(tie(it1, it2) = in_edges(x); it1 != it2; it1++)
		{
			int t = (*it1)->target();
			vd[t]--;
			assert(vd[t] >= 0);
			if(vd[t] == 0) q.push_back(t);
		}
	}
	return v;
}

vector<int> directed_graph::topological_sort()
{
	vector<int> v;

	vector<int> q;
	vector<int> vd;
	for(int i = 0; i < num_vertices(); i++)
	{
		int d = in_degree(i);
		vd.push_back(d);
		if(d == 0) q.push_back(i);
	}

	int k = 0;

	while(k < q.size())
	{
		int x = q[k++];
		v.push_back(x);

		edge_iterator it1, it2;
		for(tie(it1, it2) = out_edges(x); it1 != it2; it1++)
		{
			int t = (*it1)->target();
			vd[t]--;
			assert(vd[t] >= 0);
			if(vd[t] == 0) q.push_back(t);
		}
	}
	return v;
}

vector<int> directed_graph::topological_sort0()
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

int directed_graph::compute_in_partner(int x)
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

int directed_graph::compute_out_partner(int x)
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

int directed_graph::compute_out_equivalent_vertex(int x)
{
	// TODO: Assume that [0, gr.num_vertices() - 1] is a TPO
	vector<int> q;
	set<int> s;
	int k = 0;
	s.insert(x);
	q.push_back(x);

	int y = -1;
	while(k < q.size())
	{
		int a = q[k++];
		edge_iterator it1, it2;
		for(tie(it1, it2) = out_edges(a); it1 != it2; it1++)
		{
			int b = (*it1)->target();
			if(s.find(b) != s.end()) continue;
			s.insert(b);
			if(b > y) 
			{
				if(y != -1) q.push_back(y);
				y = b;
			}
			else
			{
				q.push_back(b);
			}
		}
	}

	if(y == -1) return -1;

	edge_iterator it1, it2;
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int a = (*it);
		if(a == x) continue;
		for(tie(it1, it2) = in_edges(a); it1 != it2; it1++)
		{
			int b = (*it1)->source();
			if(s.find(b) == s.end()) return -1;
		}
	}

	return y;
}

int directed_graph::compute_in_equivalent_vertex(int x)
{
	// TODO: Assume that [0, gr.num_vertices() - 1] is a TPO
	vector<int> q;
	set<int> s;
	int k = 0;
	s.insert(x);
	q.push_back(x);

	int y = -1;
	while(k < q.size())
	{
		int a = q[k++];
		edge_iterator it1, it2;
		for(tie(it1, it2) = in_edges(a); it1 != it2; it1++)
		{
			int b = (*it1)->source();
			if(s.find(b) != s.end()) continue;
			s.insert(b);
			if(b < y || y == -1) 
			{
				if(y != -1) q.push_back(y);
				y = b;
			}
			else
			{
				q.push_back(b);
			}
		}
	}

	if(y == -1) return -1;

	edge_iterator it1, it2;
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int a = (*it);
		if(a == x) continue;
		for(tie(it1, it2) = out_edges(a); it1 != it2; it1++)
		{
			int b = (*it1)->target();
			if(s.find(b) == s.end()) return -1;
		}
	}

	return y;
}

int directed_graph::check_nest(int x, int y, const vector<int> &tpo)
{
	set<edge_descriptor> se;
	return check_nest(x, y, se, tpo);
}

int directed_graph::check_nest(int x, int y, set<edge_descriptor> &se)
{
	vector<int> v = topological_sort();
	vector<int> tpo;
	tpo.assign(num_vertices(), -1);
	for(int i = 0; i < tpo.size(); i++)
	{
		tpo[v[i]] = i;
	}
	return check_nest(x, y, se, tpo);
}

int directed_graph::check_nest(int x, int y, set<edge_descriptor> &se, const vector<int> &tpo)
{
	vector<int> rv;
	bfs_reverse(y, rv);
	set<int> srv(rv.begin(), rv.end());

	/*
	vector<int> v = topological_sort();
	vector<int> order;
	order.assign(num_vertices(), -1);
	for(int i = 0; i < tpo.size(); i++)
	{
		order[tpo[i]] = i;
	}
	*/

	se.clear();
	set<int> sv;
	edge_iterator it1, it2;
	for(tie(it1, it2) = out_edges(x); it1 != it2; it1++)
	{
		assert((*it1)->source() == x);
		int t = (*it1)->target();
		//printf("check %x, out edge target %d, rv = %d\n", x, t, rv[t]);
		if(srv.find(t) == srv.end()) continue;
		if(sv.find(t) == sv.end()) sv.insert(t);
		if(se.find(*it1) == se.end()) se.insert(*it1);
	}

	//printf(" check (%d, %d), %lu internal vertices\n", x, y, sv.size());

	int cnt = 0;
	while(sv.size() >= 2)
	{
		int k = -1;
		for(set<int>::iterator it = sv.begin(); it != sv.end(); it++)
		{
			if(k == -1 || tpo[*it] < tpo[k])
			{
				k = *it;
			}
		}
		assert(k != -1);

		cnt++;

		for(tie(it1, it2) = in_edges(k); it1 != it2; it1++)
		{
			if(se.find(*it1) == se.end()) return -1;
		}

		for(tie(it1, it2) = out_edges(k); it1 != it2; it1++)
		{
			assert(se.find(*it1) == se.end());
			int t = (*it1)->target();
			se.insert(*it1);
			if(sv.find(t) == sv.end()) sv.insert(t);
		}
		sv.erase(k);
	}

	if(sv.size() == 0) return -1;
	assert(sv.size() == 1);
	if(*(sv.begin()) == y) return cnt;
	else return -1;
}

int directed_graph::draw(const string &file, const MIS &mis, const MES &mes, double len)
{
	vector<int> tp;
	for(int i = 0; i < num_vertices(); i++) tp.push_back(i);
	return draw(file, mis, mes, len, tp);
}

int directed_graph::draw(const string &file, const MIS &mis, const MES &mes, double len, const vector<int> &tp)
{
	ofstream fout(file.c_str());
	if(fout.fail())
	{
		printf("open file %s error.\n", file.c_str());
		return 0;
	}

	draw_header(fout);

	fout<<"\\def\\len{"<<len<<"cm}\n";

	// draw file name
	fout<<"\\node[draw, thick, red] at (1.6 * \\len, 0.58 * \\len) {"<<file.c_str()<<"};\n";

	// draw vertices
	char sx[1024];
	char sy[1024];
	double pos = 0;

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
			//if(i + 1 == j) bend = 0;

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

