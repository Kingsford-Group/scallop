/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "undirected_graph.h"
#include "draw.h"

#include <climits>
#include <cstdio>
#include <cassert>
#include <tuple>
#include <algorithm>

using namespace std;

undirected_graph::undirected_graph()
{}

undirected_graph::undirected_graph(const undirected_graph &gr)
{
	copy(gr);
}

undirected_graph& undirected_graph::operator=(const undirected_graph &gr)
{
	copy(gr);
	return (*this);
}

undirected_graph::~undirected_graph()
{
}

edge_descriptor undirected_graph::add_edge(int s, int t)
{
	assert(s >= 0 && s < vv.size());
	assert(t >= 0 && t < vv.size());
	edge_base *e = new edge_base(s, t);
	assert(se.find(e) == se.end());
	se.insert(e);
	vv[s]->add_out_edge(e);
	vv[t]->add_out_edge(e);
	return e;
}

int undirected_graph::remove_edge(edge_descriptor e)
{
	if(se.find(e) == se.end()) return -1;
	vv[e->source()]->remove_out_edge(e);
	vv[e->target()]->remove_out_edge(e);
	delete e;
	se.erase(e);
	return 0;
}

int undirected_graph::remove_edge(int s, int t)
{
	vector<edge_base*> v;
	PEEI p = vv[s]->out_edges();
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		int ss = (*it)->source();
		int tt = (*it)->target();
		assert(ss == s || tt == s);
		if(ss == s && tt != t) continue;
		if(tt == s && ss != t) continue;
		v.push_back(*it);
	}
	for(int i = 0; i < v.size(); i++)
	{
		remove_edge(v[i]);
	}
	return 0;
}

PEEI undirected_graph::out_edges(int x)
{
	PEEI p = vv[x]->out_edges();
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		if((*it)->source() != x) (*it)->swap();
	}
	return p;
}

vector<int> undirected_graph::assign_connected_components()
{
	vector<bool> mm;
	vector<int> vv;
	mm.assign(num_vertices(), false);
	vv.assign(num_vertices(), -1);

	int cc = 0;
	for(int i = 0; i < num_vertices(); i++)
	{
		if(mm[i] == true) continue;
		vector<int> v;
		bfs(i, v);
		for(int k = 0; k < v.size(); k++) mm[v[k]] = true;
		for(int k = 0; k < v.size(); k++) vv[v[k]] = cc;
		cc++;
	}
	return vv;
}

vector< set<int> > undirected_graph::compute_connected_components()
{
	vector<bool> m;
	m.assign(num_vertices(), false);
	vector< set<int> > vv;

	for(int i = 0; i < num_vertices(); i++)
	{
		if(m[i] == true) continue;
		vector<int> v;
		bfs(i, v);
		set<int> s(v.begin(), v.end());
		assert(s.find(i) != s.end());
		vv.push_back(s);
		for(int k = 0; k < v.size(); k++) m[v[k]] = true;
	}
	return vv;
}

bool undirected_graph::intersect(edge_descriptor ex, edge_descriptor ey)
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

int undirected_graph::draw(const string &file, const MIS &mis, const MES &mes, double len)
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
			if(i >= j) continue;
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
			fout<<"\\draw["<< line.c_str() <<"\\colx, bend right = "<< bend <<"] ("<<sx<<") to node[draw, fill=white] {";
			fout<< s.c_str() <<"} ("<<sy<<");\n";
		}
	}

	draw_footer(fout);

	fout.close();
	return 0;
}

