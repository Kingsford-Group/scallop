#include "nested_graph.h"
#include <algorithm>

nested_graph::nested_graph()
{
}

nested_graph::~nested_graph()
{}

int nested_graph::solve(directed_graph &gr)
{
	init(gr);
	build_nests(gr);
	get_edge_indices(i2e, e2i);
	build_partial_order();
	build_partners(gr);
	build_parents();
	build_docking();
	test_linking(gr);
	return 0;
}

int nested_graph::init(directed_graph &gr)
{
	clear();
	for(int i = 0; i < gr.num_vertices(); i++) 
	{
		add_vertex();
	}
	return 0;
}

int nested_graph::build_nests(directed_graph &gr)
{
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		if(gr.degree(i) == 0) continue;
		for(int j = 0; j < gr.num_vertices(); j++)
		{
			if(gr.degree(j) == 0) continue;
			int b = gr.check_nest(i, j);
			if(b == -1) continue;
			add_edge(i, j);
		}
	}
	assert(check_nested() == true);
	return 0;
}

int nested_graph::build_partners(directed_graph &gr)
{
	partners.clear();
	partners.assign(gr.num_vertices(), PI(-1, -1));
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		if(gr.degree(i) == 0) continue;
		int ip = gr.compute_in_partner(i);
		int op = gr.compute_out_partner(i);
		partners[i] = PI(ip, op);
	}
	return 0;
}

int nested_graph::build_partial_order()
{
	order.assign(num_vertices(), -1);
	vector<int> v = topological_sort();
	for(int i = 0; i < v.size(); i++) 
	{
		order[v[i]] = i;
	}
	return 0;
}

int nested_graph::build_parents()
{
	parents.assign(num_edges(), -1);
	vector<int> v = topological_sort();
	for(int i = 0; i < v.size(); i++)
	{
		int x = v[i];
		if(degree(x) == 0) continue;
		build_parents(x);
	}
	return 0;
}

int nested_graph::build_parents(int x)
{
	edge_iterator it1, it2;
	int k = x;
	int e = -1;
	for(tie(it1, it2) = in_edges(x); it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		assert(t == x);
		if(order[s] < order[k])
		{
			k = s;
			e = e2i[*it1];
		}
	}
	int p = -1;
	if(e != -1) p = parents[e];

	vector<PI> pp;
	for(tie(it1, it2) = out_edges(x); it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		assert(s == x);
		pp.push_back(PI(order[t], e2i[*it1]));
	}

	sort(pp.begin(), pp.end());

	for(int i = pp.size() - 1; i >= 0; i--)
	{
		int q = p;
		int ei = pp[i].second;
		int ii = i2e[ei]->target();
		for(int j = i + 1; j < pp.size(); j++)
		{
			int ej = pp[j].second;
			int jj = i2e[ej]->target();
			if(check_path(ii, jj) == false) continue;
			q = ej;
			break;
		}
		parents[ei] = q;
	}
	return 0;
}

int nested_graph::compute_parent(edge_descriptor e)
{
	int ss = e->source();
	int tt = e->target();

	vector<PI> pp;
	edge_iterator it1, it2;
	for(tie(it1, it2) = out_edges(ss); it1 != it2; it1++)
	{
		if((*it1) == e) continue;
		int s = (*it1)->source();
		int t = (*it1)->target();
		assert(s == ss);
		pp.push_back(PI(order[t], e2i[*it1]));
	}

	sort(pp.begin(), pp.end());

	for(int i = 0; i < pp.size(); i++)
	{
		int ei = pp[i].second;
		int t = i2e[ei]->target();
		if(t == tt) continue;
		if(order[t] <= order[tt]) continue;
		if(check_path(tt, t) == true)
		{
			return ei;
		}
	}

	int k = ss;
	int ee = -1;
	for(tie(it1, it2) = in_edges(ss); it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		assert(t == ss);
		if(order[s] < order[k])
		{
			k = s;
			ee = e2i[*it1];
		}
	}
	int p = -1;
	if(ee != -1) p = parents[ee];

	return p;
}

int nested_graph::build_docking()
{
	docking.assign(num_edges(), false);

	vector<int> v;
	v = bfs_docking_forward(0, -1);
	for(int k = 0; k < v.size(); k++) docking[v[k]] = true;

	v = bfs_docking_backward(num_vertices() - 1, -1);
	for(int k = 0; k < v.size(); k++) docking[v[k]] = true;

	for(int i = 0; i < num_edges(); i++)
	{
		int s = i2e[i]->source();
		int t = i2e[i]->target();
		v = bfs_docking_forward(s, i);
		for(int k = 0; k < v.size(); k++) docking[v[k]] = true;
		v = bfs_docking_backward(t, i);
		for(int k = 0; k < v.size(); k++) docking[v[k]] = true;
	}
	return 0;
}

vector<int> nested_graph::bfs_docking_forward(int s, int pp)
{
	set<int> sv;
	vector<int> v;
	edge_iterator it1, it2;
	for(tie(it1, it2) = out_edges(s); it1 != it2; it1++)
	{
		int e = e2i[*it1];
		if(parents[e] != pp) continue;
		//docking[e] = true;
		v.push_back(e);
		int t = (*it1)->target();
		assert(sv.find(t) == sv.end());
		sv.insert(t);
	}

	int p = 0;
	while(p < v.size())
	{
		int ee = v[p];
		p++;

		int ss = i2e[ee]->source();
		int tt = i2e[ee]->target();

		if(partners[tt].first == -1 || partners[tt].second == -1) continue;
		if(partners[tt].first != ss) continue;

		for(tie(it1, it2) = out_edges(tt); it1 != it2; it1++)
		{
			int e = e2i[*it1];
			if(parents[e] != pp) continue;
			if(partners[tt].second != (*it1)->target()) continue;
			//assert(docking[e] == false);
			//docking[e] = true;
			v.push_back(e);
			int t = (*it1)->target();
			//assert(sv.find(t) == sv.end());
			sv.insert(s);
		}
	}
	return v;
}

vector<int> nested_graph::bfs_docking_backward(int t, int pp)
{
	set<int> sv;
	vector<int> v;
	edge_iterator it1, it2;
	for(tie(it1, it2) = in_edges(t); it1 != it2; it1++)
	{
		int e = e2i[*it1];
		if(parents[e] != pp) continue;
		//if(docking[e] == true) continue;
		//docking[e] = true;
		v.push_back(e);
		int s = (*it1)->source();
		assert(sv.find(s) == sv.end());
		sv.insert(s);
	}

	int p = 0;
	while(p < v.size())
	{
		int ee = v[p];
		p++;

		int ss = i2e[ee]->source();
		int tt = i2e[ee]->target();

		if(partners[ss].first == -1 || partners[ss].second == -1) continue;
		if(partners[ss].second != tt) continue;

		for(tie(it1, it2) = in_edges(ss); it1 != it2; it1++)
		{
			int e = e2i[*it1];
			if(parents[e] != pp) continue;
			if(partners[ss].first != (*it1)->source()) continue;
			//if(docking[e] == true) continue;
			//docking[e] = true;
			v.push_back(e);
			int s = (*it1)->source();
			//assert(sv.find(s) == sv.end());
			sv.insert(s);
		}
	}
	return v;
}

int nested_graph::compute_lca(int x, int y, vector<int> &xx, vector<int> &yy)
{
	vector<int> xv;
	vector<int> yv;

	xv.push_back(x);
	yv.push_back(y);
	int p = x;
	while(p != -1)
	{
		p = parents[p];
		xv.push_back(p);
	}
	p = y;
	while(p != -1)
	{
		p = parents[p];
		yv.push_back(p);
	}

	xx.clear();
	yy.clear();
	for(int i = 0; i < xv.size(); i++)
	{
		for(int j = 0; j < yv.size(); j++)
		{
			if(xv[i] != yv[j]) continue;
			xx.assign(xv.begin(), xv.begin() + i + 1);
			yy.assign(yv.begin(), yv.begin() + j + 1);
			return xv[i];
		}
	}
	assert(false);
}

bool nested_graph::dock(int e, int p)
{
	int x = e;
	while(x != p)
	{
		if(docking[x] == false) return false;
		x = parents[x];
	}
	return true;
}

PEB nested_graph::add_extra_edge(int s, int t)
{
	PEB p = edge(s, t);
	if(p.second == true) return p;

	edge_descriptor e = add_edge(s, t);
	int n = i2e.size();
	i2e.push_back(e);
	e2i.insert(PEI(e, n));

	int pp = compute_parent(e);
	parents.push_back(pp);

	docking.push_back(false);
	if(i2e[pp]->source() == s) docking[n] = true;
	if(i2e[pp]->target() == t) docking[n] = true;

	return PEB(e, false);
}

int nested_graph::remove_extra_edge(PEB p)
{
	if(p.second == true) return 0;
	int i = e2i[p.first];
	assert(i == e2i.size() - 1);
	remove_edge(p.first);

	e2i.erase(p.first);
	i2e.pop_back();

	parents.pop_back();
	docking.pop_back();

	return 0;
}

bool nested_graph::link(int xs, int xt, int ys, int yt, vector<PI> &p)
{
	p.clear();
	if(xs == yt) return true;
	if(xt == ys) return true;
	if(xs == ys) return false;
	if(xt == yt) return false;

	PEB xb = add_extra_edge(xs, xt);
	PEB yb = add_extra_edge(ys, yt);

	edge_descriptor xe = xb.first;
	edge_descriptor ye = yb.first;
	int xi = e2i[xe];
	int yi = e2i[ye];

	if(intersect(xe, ye) == true)
	{
		remove_extra_edge(yb);
		remove_extra_edge(xb);
		return false;
	}

	assert(check_nested() == true);

	vector<int> xv;
	vector<int> yv;
	int lca = compute_lca(xi, yi, xv, yv);

	if(xv.size() < 2 || yv.size() < 2)
	{
		remove_extra_edge(yb);
		remove_extra_edge(xb);
		return false;
	}

	int xu = xv[xv.size() - 2];
	int yu = yv[yv.size() - 2];

	bool xd = dock(xi, xu);
	bool yd = dock(yi, yu);

	/*
	printf("edges %d:(%d, %d) and %d:(%d, %d), parents = (%d, %d), lca = %d, sub-lca = (%d, %d), dock = (%c, %c)\n",
			xi, xs, xt, yi, ys, yt, parents[xi], parents[yi], lca, xu, yu, xd ? 'T' : 'F', yd ? 'T' : 'F'); 
	*/

	if(xd == false || yd == false)
	{
		remove_extra_edge(yb);
		remove_extra_edge(xb);
		return false;
	}

	assert(xu >= 0 && xu < num_edges());
	assert(yu >= 0 && yu < num_edges());

	vector<int> xx = bfs_docking_forward(i2e[xu]->target(), lca);
	vector<int> yy = bfs_docking_forward(i2e[yu]->target(), lca);

	for(int i = 0; i < xx.size(); i++)
	{
		if(xx[i] != yu) continue;
		remove_extra_edge(yb);
		remove_extra_edge(xb);
		return true;
	}

	for(int i = 0; i < yy.size(); i++)
	{
		if(yy[i] != xu) continue;
		remove_extra_edge(yb);
		remove_extra_edge(xb);
		return true;
	}

	remove_extra_edge(yb);
	remove_extra_edge(xb);
	return false;
}

int nested_graph::test_linking(directed_graph &gr)
{
	edge_iterator i1, i2;
	edge_iterator j1, j2;
	for(tie(i1, i2) = gr.edges(); i1 != i2; i1++)
	{
		int xs = (*i1)->source();
		int xt = (*i1)->target();
		for(tie(j1, j2) = gr.edges(); j1 != j2; j1++)
		{
			int ys = (*j1)->source();
			int yt = (*j1)->target();
			vector<PI> p;
			bool b = link(xs, xt, ys, yt, p);

			printf("link (%d, %d) and (%d, %d) = %c\n", xs, xt, ys, yt, b ? 'T' : 'F');
		}
	}
	return 0;
}

int nested_graph::draw(const string &file) 
{
	MIS mis;
	MES mes;
	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(); it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		string ss;
		char buf[1024];
		sprintf(buf, "%d(", e2i[*it1]);
		ss = string(buf);

		bool b = false;
		if(partners[s].second == t && partners[s].first >= 0)
		{
			sprintf(buf, "%d", s);
			ss += string(buf);
			b = true;
		}
		if(partners[t].first == s && partners[t].second >= 0)
		{
			sprintf(buf, "%d", t);
			if(b == false) ss += string(buf);
			else ss += (string(",") + string(buf));
		}
		ss += string(")");

		sprintf(buf, "%d", parents[e2i[*it1]]);
		ss += string(buf);
		
		if(docking[e2i[*it1]] == true) ss += "T";
		else ss += "F";

		mes.insert(PES(*it1, ss));
	}
	directed_graph::draw(file, mis, mes, 3.5);
	return 0;
}

int nested_graph::print()
{
	// print parents:
	for(int i = 0; i < parents.size(); i++)
	{
		int s = i2e[i]->source();
		int t = i2e[i]->target();
		printf("edge (%d, %d), index = %d, parent = %d, docking = %c\n", s, t, i, parents[i], docking[i] ? 'T' : 'F');
	}
	return 0;
}
