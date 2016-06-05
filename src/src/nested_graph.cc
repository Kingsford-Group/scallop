#include "nested_graph.h"
#include <algorithm>

nested_graph::nested_graph()
{
}

nested_graph::~nested_graph()
{}

int nested_graph::build(directed_graph &gr)
{
	clear();
	init(gr);
	build_nests(gr);
	get_edge_indices(i2e, e2i);
	build_partial_order();
	build_parents();
	build_partners(gr);
	//test_linking(gr);
	return 0;
}

int nested_graph::clear()
{
	directed_graph::clear();
	i2e.clear();
	e2i.clear();
	order.clear();
	partners.clear();
	parents.clear();
	return 0;
}

int nested_graph::init(directed_graph &gr)
{
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

		// checking consistent with left- right-siblings
		if(ip < 0 || op < 0) continue;
		PEB p1 = edge(ip, i);
		PEB p2 = edge(i, op);
		assert(p1.second && p2.second);
		assert(e2i[p1.first] == get_left_sibling(e2i[p2.first]));
		assert(e2i[p2.first] == get_right_sibling(e2i[p1.first]));
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

int nested_graph::get_left_sibling(int e)
{
	int p = parents[e];
	int s = i2e[e]->source();
	edge_iterator it1, it2;
	for(tie(it1, it2) = in_edges(s); it1 != it2; it1++)
	{
		int ee = e2i[*it1];
		int pp = parents[ee];
		if(pp == p) return ee;
	}
	return -1;
}

int nested_graph::get_right_sibling(int e)
{
	int p = parents[e];
	int t = i2e[e]->target();
	edge_iterator it1, it2;
	for(tie(it1, it2) = out_edges(t); it1 != it2; it1++)
	{
		int ee = e2i[*it1];
		int pp = parents[ee];
		if(pp == p) return ee;
	}
	return -1;
}

bool nested_graph::bfs_docking_forward(int e0, int t0, int p0, vector<int> &r)
{
	r.clear();
	int s = i2e[e0]->source();
	int t = i2e[e0]->target();
	// TODO, while is not necessary
	while(true)
	{
		if(t == t0) return true;
		if(partners[t].first == -1 || partners[t].second == -1) return false;
		if(partners[t].first != s) return false;
		bool b = false;
		edge_iterator it1, it2;
		for(tie(it1, it2) = out_edges(t); it1 != it2; it1++)
		{
			int ee = e2i[*it1];
			int tt = (*it1)->target();
			if(parents[ee] != p0) continue;
			if(partners[t].second != tt) continue;
			b = true;
			r.push_back(t);
			s = t;
			t = tt;
		}
		if(b == false) return false;
	}
	return true;
}

bool nested_graph::bfs_docking_backward(int e0, int s0, int p0, vector<int> &r)
{
	r.clear();
	int s = i2e[e0]->source();
	int t = i2e[e0]->target();
	while(true)
	{
		if(s == s0) return true;
		if(partners[s].first == -1 || partners[s].second == -1) return false;
		if(partners[s].second != t) return false;
		bool b = false;
		edge_iterator it1, it2;
		for(tie(it1, it2) = in_edges(s); it1 != it2; it1++)
		{
			int ee = e2i[*it1];
			int ss = (*it1)->source();
			if(parents[ee] != p0) continue;
			if(partners[s].first != ss) continue;
			b = true;
			r.push_back(s);
			t = s;
			s = ss;
		}
		if(b == false) return false;
	}
	return true;
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

int nested_graph::dock(int e, int p, vector<PI> &v)
{
	v.clear();
	int x = e;
	int d = 0;
	while(x != p)
	{
		int s = i2e[x]->source();
		int t = i2e[x]->target();

		int xx = parents[x];
		int ss = i2e[xx]->source();
		int tt = i2e[xx]->target();

		vector<int> vt;
		vector<int> vs;
		bool bt = bfs_docking_forward(x, tt, xx, vt);
		bool bs = bfs_docking_backward(x, ss, xx, vs);

		if(bt == true)
		{
			if(d < 0) v.push_back(PI(s, t));
			for(int k = 0; k < vt.size(); k++) v.push_back(PI(vt[k], -1));
			d = 1;
		}
		else if(bs == true)
		{
			if(d > 0) v.push_back(PI(s, t));
			for(int k = 0; k < vs.size(); k++) v.push_back(PI(vs[k], -1));
			d = -1;
		}
		else return -2;

		x = xx;
	}
	return d;
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

	return 0;
}

bool nested_graph::link(int xs, int xt, int ys, int yt, vector<PI> &xp, vector<PI> &yp)
{
	xp.clear();
	yp.clear();

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

	bool b1 = check_path(i2e[xu], i2e[yu]);
	bool b2 = check_path(i2e[yu], i2e[xu]);

	assert(b1 || b2);

	assert(xu >= 0 && xu < num_edges());
	assert(yu >= 0 && yu < num_edges());

	int xd = dock(xi, xu, xp);
	int yd = dock(yi, yu, yp);

	/*
	printf("edges %d:(%d, %d) and %d:(%d, %d), parents = (%d, %d), lca = %d, sub-lca = (%d, %d), dock = (%c, %c)\n",
			xi, xs, xt, yi, ys, yt, parents[xi], parents[yi], lca, xu, yu, xd ? 'T' : 'F', yd ? 'T' : 'F'); 
	*/

	if(xd == -2 || yd == -2)
	{
		remove_extra_edge(yb);
		remove_extra_edge(xb);
		return false;
	}

	int xus = i2e[xu]->source();
	int xut = i2e[xu]->target();
	int yus = i2e[yu]->source();
	int yut = i2e[yu]->target();

	vector<int> vxyf, vxyb, vyxf, vyxb;
	bool xyf = bfs_docking_forward(xu, yus, lca, vxyf);
	bool yxf = bfs_docking_forward(yu, xus, lca, vyxf);
	bool xyb = bfs_docking_backward(xu, yut, lca, vxyb);
	bool yxb = bfs_docking_backward(yu, xut, lca, vyxb);


	if(xyf == true)
	{
		if(xd == -1) xp.push_back(PI(xus, xut));
		if(yd == +1) yp.push_back(PI(yus, yut));
		for(int k = 0; k < vxyf.size(); k++) xp.push_back(PI(vxyf[k], -1));
		remove_extra_edge(yb);
		remove_extra_edge(xb);
		return true;
	}
	else if(yxf == true)
	{
		if(yd == -1) yp.push_back(PI(yus, yut));
		if(xd == +1) xp.push_back(PI(xus, xut));
		for(int k = 0; k < vyxf.size(); k++) yp.push_back(PI(vyxf[k], -1));
		remove_extra_edge(yb);
		remove_extra_edge(xb);
		return true;
	}
	else if(xyb == true)
	{
		if(xd == -1) xp.push_back(PI(xus, xut));
		if(yd == +1) yp.push_back(PI(yus, yut));
		for(int k = 0; k < vxyb.size(); k++) yp.push_back(PI(vxyb[k], -1));
		remove_extra_edge(yb);
		remove_extra_edge(xb);
		return true;
	}
	else if(yxb == true)
	{
		if(yd == -1) yp.push_back(PI(yus, yut));
		if(xd == +1) xp.push_back(PI(xus, xut));
		for(int k = 0; k < vyxb.size(); k++) xp.push_back(PI(vyxb[k], -1));
		remove_extra_edge(yb);
		remove_extra_edge(xb);
		return true;
	}
	else
	{
		remove_extra_edge(yb);
		remove_extra_edge(xb);
		return false;
	}
}

int nested_graph::test_linking(directed_graph &gr)
{
	edge_iterator i1, i2;
	edge_iterator j1, j2;
	set<PI> s1;
	for(tie(i1, i2) = gr.edges(); i1 != i2; i1++)
	{
		int xs = (*i1)->source();
		int xt = (*i1)->target();
		if(s1.find(PI(xs, xt)) != s1.end()) continue;
		s1.insert(PI(xs, xt));
		set<PI> s2;
		for(tie(j1, j2) = gr.edges(); j1 != j2; j1++)
		{
			int ys = (*j1)->source();
			int yt = (*j1)->target();
			if(s2.find(PI(ys, yt)) != s2.end()) continue;
			s2.insert(PI(ys, yt));
			vector<PI> xp, yp;
			bool b = link(xs, xt, ys, yt, xp, yp);

			printf("link (%d, %d) and (%d, %d) = %c\n", xs, xt, ys, yt, b ? 'T' : 'F');
			printf(" X-action:");
			for(int k = 0; k < xp.size(); k++) printf(" (%d, %d), ", xp[k].first, xp[k].second);
			printf("\n Y-action:");
			for(int k = 0; k < yp.size(); k++) printf(" (%d, %d), ", yp[k].first, yp[k].second);
			printf("\n");
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
		printf("edge (%d, %d), index = %d, parent = %d\n", s, t, i, parents[i]);
	}
	return 0;
}
