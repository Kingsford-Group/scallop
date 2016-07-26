#include "nested_graph.h"
#include "util.h"
#include <algorithm>

nested_graph::nested_graph()
{
}

nested_graph::nested_graph(directed_graph &gr)
{
	build(gr);
}

nested_graph::~nested_graph()
{}

int nested_graph::build(directed_graph &gr)
{
	clear();
	init(gr);
	build_nests0(gr);
	build_partial_order();
	build_parents();
	build_partners(gr);
	return 0;
}

int nested_graph::clear()
{
	directed_graph::clear();
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
	vector< set<int> > vs;
	vector< set<int> > vt;
	vs.resize(gr.num_vertices());
	vt.resize(gr.num_vertices());
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		if(gr.degree(i) == 0) continue;
		vector<int> v;
		set<int> s;
		gr.bfs(i, v);
		for(int k = 0; k < v.size(); k++)
		{
			s.insert(v[k]);
		}
		vs[i] = s;

		v.clear();
		s.clear();
		gr.bfs_reverse(i, v);
		for(int k = 0; k < v.size(); k++)
		{
			s.insert(v[k]);
		}
		vt[i] = s;
	}

	vector<int> v = gr.topological_sort();
	assert(v.size() == gr.num_vertices());
	vector<int> tpo;
	tpo.assign(gr.num_vertices(), -1);
	for(int i = 0; i < v.size(); i++)
	{
		tpo[v[i]] = i;
	}

	for(int i = 0; i < gr.num_vertices(); i++)
	{
		if(gr.degree(i) == 0) continue;
		for(int j = 0; j < gr.num_vertices(); j++)
		{
			if(gr.degree(j) == 0) continue;

			// for verifying
			//int f = gr.check_nest(i, j);
			//if(f == -1) continue;

			bool b = verify_minimal_nest(gr, tpo, vs, vt, i, j);
			if(b == false) continue;

			add_edge(i, j);
		}
	}
	//assert(check_nested() == true);
	return 0;
}

bool nested_graph::verify_minimal_nest(directed_graph &gr, const vector<int> &tpo, vector< set<int> >&vs, vector< set<int> > &vt, int i , int j)
{
	vector<int> vv;
	vv.resize(gr.num_vertices());

	vector<int>::iterator iv = set_intersection(vs[i].begin(), vs[i].end(), vt[j].begin(), vt[j].end(), vv.begin());
	set<int> s(vv.begin(), iv);

	/*
	printf("vertices = ");
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		printf("%d:%d ", *it, tpo[*it]);
	}
	printf("\n");
	*/
	if(s.size() <= 1) return false;

	// verify that it is a nest
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int x = *it;
		if(x == i || x == j) continue;
		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.in_edges(x); it1 != it2; it1++)
		{
			int y = (*it1)->source();
			if(s.find(y) == s.end()) return false;
		}
		for(tie(it1, it2) = gr.out_edges(x); it1 != it2; it1++)
		{
			int y = (*it1)->target();
			if(s.find(y) == s.end()) return false;
		}
	}

	// verify that it is a minimal nest
	set<int> ss;
	vector<int> q;
	int k = 0;
	ss.insert(i);
	q.push_back(i);
	int m = -1;
	while(true)
	{
		if(k == q.size()) 
		{
			if(m == j) return true;
			else return false;
		}
		int x = q[k++];
		
		//printf("check x = %d, m = %d, k = %d, q = %lu, ss = %lu\n", x, m, k, q.size(), ss.size()); 
		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.out_edges(x); it1 != it2; it1++)
		{
			int y = (*it1)->target();
			//printf(" adjacent to %d\n", y);
			if(s.find(y) == s.end()) continue;
			if(ss.find(y) != ss.end()) continue;
			if(y == m) continue;
			ss.insert(x);
			if(m == -1 || tpo[y] > tpo[m])
			{
				if(m != -1) q.push_back(m);
				m = y;
			}
			else
			{
				q.push_back(y);
			}
		}
	}
	assert(false);
}

int nested_graph::build_nests0(directed_graph &gr)
{
	vector<int> v1;
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		if(gr.degree(i) == 0) continue;
		v1.push_back(i);
	}

	vector<int> v = gr.topological_sort();
	vector<int> tpo;
	tpo.assign(num_vertices(), -1);
	for(int i = 0; i < tpo.size(); i++)
	{
		tpo[v[i]] = i;
	}

	int c = 0;
	for(int i = 0; i < v1.size(); i++)
	{
		for(int j = 0; j < v1.size(); j++)
		{
			int b = gr.check_nest(v1[i], v1[j], tpo);
			if(b == -1) continue;
			add_edge(v1[i], v1[j]);
			c++;
		}
	}
	//printf("building nested graph takes %d check_nest\n", c);
	//assert(check_nested() == true);
	return 0;
}

int nested_graph::build_partners(directed_graph &gr)
{
	partners.clear();
	partners.assign(gr.num_vertices(), PEE(null_edge, null_edge));
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		if(gr.degree(i) == 0) continue;
		int ip = gr.compute_in_partner(i);
		int op = gr.compute_out_partner(i);
		if(ip >= 0)
		{
			PEB p1 = edge(ip, i);
			assert(p1.second == true);
			assert(ip == p1.first->source());
			partners[i].first = p1.first;
		}
		if(op >= 0)
		{
			PEB p2 = edge(i, op);
			assert(p2.second == true);
			assert(op == p2.first->target());
			partners[i].second = p2.first;
		}
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
	parents.clear();
	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(); it1 != it2; it1++)
	{
		parents.insert(PEE(*it1, null_edge));
	}

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
	edge_descriptor e = null_edge;
	for(tie(it1, it2) = in_edges(x); it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		assert(t == x);
		if(order[s] < order[k])
		{
			k = s;
			e = *it1;
		}
	}
	edge_descriptor p = null_edge;
	if(e != null_edge) p = parents[e];

	vector<PIE> pp;
	for(tie(it1, it2) = out_edges(x); it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		assert(s == x);
		pp.push_back(PIE(order[t], *it1));
	}

	sort(pp.begin(), pp.end());

	for(int i = pp.size() - 1; i >= 0; i--)
	{
		edge_descriptor q = p;
		edge_descriptor ei = pp[i].second;
		int ii = ei->target();
		for(int j = i + 1; j < pp.size(); j++)
		{
			edge_descriptor ej = pp[j].second;
			int jj = ej->target();
			if(check_path(ii, jj) == false) continue;
			q = ej;
			break;
		}
		parents[ei] = q;
	}
	return 0;
}

edge_descriptor nested_graph::compute_parent(edge_descriptor e)
{
	int ss = e->source();
	int tt = e->target();

	vector<PIE> pp;
	edge_iterator it1, it2;
	for(tie(it1, it2) = out_edges(ss); it1 != it2; it1++)
	{
		if((*it1) == e) continue;
		int s = (*it1)->source();
		int t = (*it1)->target();
		assert(s == ss);
		pp.push_back(PIE(order[t], *it1));
	}

	sort(pp.begin(), pp.end());

	for(int i = 0; i < pp.size(); i++)
	{
		edge_descriptor ei = pp[i].second;
		int t = ei->target();
		if(t == tt) continue;
		if(order[t] <= order[tt]) continue;
		if(check_path(tt, t) == true) return ei;
	}

	int k = ss;
	edge_descriptor ee = null_edge;
	for(tie(it1, it2) = in_edges(ss); it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		assert(t == ss);
		if(order[s] < order[k])
		{
			k = s;
			ee = *it1;
		}
	}
	edge_descriptor p = null_edge;
	if(ee != null_edge) p = parents[ee];

	return p;
}

edge_descriptor nested_graph::get_left_sibling(edge_descriptor e)
{
	edge_descriptor p = parents[e];
	int s = e->source();
	edge_iterator it1, it2;
	for(tie(it1, it2) = in_edges(s); it1 != it2; it1++)
	{
		if(parents[*it1] == p) return (*it1);
	}
	return null_edge;
}

edge_descriptor nested_graph::get_right_sibling(edge_descriptor e)
{
	edge_descriptor p = parents[e];
	int t = e->target();
	edge_iterator it1, it2;
	for(tie(it1, it2) = out_edges(t); it1 != it2; it1++)
	{
		if(parents[*it1] == p) return (*it1);
	}
	return null_edge;
}

int nested_graph::move_left(edge_descriptor e)
{
	edge_descriptor x = get_left_sibling(e);
	if(x == null_edge) return -1;
	int s = e->source();
	int t = e->target();
	if(partners[s].first != x) return -1;
	if(partners[s].second != e) return -1;
	return s;
}

int nested_graph::move_right(edge_descriptor e)
{
	edge_descriptor x = get_right_sibling(e);
	if(x == null_edge) return -1;
	int s = e->source();
	int t = e->target();
	if(partners[t].first != e) return -1;
	if(partners[t].second != x) return -1;
	return t;
}

bool nested_graph::dock_left(edge_descriptor e, vector<int> &r)
{
	r.clear();
	edge_descriptor p = parents[e];
	int d = (p == null_edge) ? 0 : p->source();

	edge_descriptor ee = e;
	while(true)
	{
		if(ee->source() == d) break;
		int m = move_left(ee);
		if(m < 0) return false;
		r.push_back(m);
		ee = get_left_sibling(ee);
	}
	return true;
}

bool nested_graph::dock_right(edge_descriptor e, vector<int> &r)
{
	r.clear();
	edge_descriptor p = parents[e];
	int d = (p == null_edge) ? (num_vertices() - 1) : p->target();

	edge_descriptor ee = e;
	while(true)
	{
		if(ee->target() == d) break;
		int m = move_right(ee);
		if(m < 0) return false;
		r.push_back(m);
		ee = get_right_sibling(ee);
	}
	return true;
}

edge_descriptor nested_graph::compute_lca(edge_descriptor x, edge_descriptor y)
{
	VE xv;
	VE yv;

	xv.push_back(x);
	yv.push_back(y);

	edge_descriptor p = x;
	while(p != null_edge)
	{
		p = parents[p];
		xv.push_back(p);
	}
	p = y;
	while(p != null_edge)
	{
		p = parents[p];
		yv.push_back(p);
	}

	for(int i = 0; i < xv.size(); i++)
	{
		for(int j = 0; j < yv.size(); j++)
		{
			if(xv[i] != yv[j]) continue;
			return xv[i];
		}
	}
	assert(false);
}

int nested_graph::dock(edge_descriptor e, edge_descriptor p, vector<PI> &v)
{
	// assume p is on the path from e to root
	v.clear();
	edge_descriptor x = e;
	int d = 0;
	while(x != p)
	{
		int s = x->source();
		int t = x->target();
		vector<int> vs;
		vector<int> vt;
		bool bs = dock_left(x, vs);
		bool bt = dock_right(x, vt);

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

		x = parents[x];
	}
	return d;
}

PEB nested_graph::add_extra_edge(int s, int t)
{
	PEB p = edge(s, t);
	if(p.second == true) return p;
	edge_descriptor e = add_edge(s, t);
	edge_descriptor x = compute_parent(e);
	assert(parents.find(e) == parents.end());
	parents.insert(PEE(e, x));
	return PEB(e, false);
}

int nested_graph::remove_extra_edge(PEB p)
{
	if(p.second == true) return 0;
	remove_edge(p.first);
	assert(parents.find(p.first) != parents.end());
	parents.erase(p.first);
	return 0;
}

bool nested_graph::link(int xs, int xt, int ys, int yt)
{
	vector<PI> xp, yp;
	return link(xs, xt, ys, yt, xp, yp);
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

	if(intersect(xe, ye) == true)
	{
		remove_extra_edge(yb);
		remove_extra_edge(xb);
		return false;
	}

	//assert(check_nested() == true);

	edge_descriptor lca = compute_lca(xe, ye);

	edge_descriptor xu = xe;
	edge_descriptor yu = ye;

	while(parents[xu] != lca) xu = parents[xu];
	while(parents[yu] != lca) yu = parents[yu];

	int xd = dock(xe, xu, xp);
	int yd = dock(ye, yu, yp);

	if(xd == -2 || yd == -2)
	{
		remove_extra_edge(yb);
		remove_extra_edge(xb);
		return false;
	}

	assert(check_path(xu, yu) == true);
	
	if(xd == -1) xp.push_back(PI(xu->source(), xu->target()));
	if(yd == +1) yp.push_back(PI(yu->source(), yu->target()));

	while(true)
	{
		if(xu->target() == yu->source()) break;
		int mx = move_right(xu);
		int my = move_left(yu);
		if(mx >= 0)
		{
			xp.push_back(PI(mx, -1));
			xu = get_right_sibling(xu);
		}
		else if(my >= 0)
		{
			yp.push_back(PI(my, -1));
			yu = get_left_sibling(yu);
		}
		else
		{
			remove_extra_edge(yb);
			remove_extra_edge(xb);
			return false;
		}
	}

	remove_extra_edge(yb);
	remove_extra_edge(xb);
	return true;
}

int nested_graph::draw(const string &file) 
{
	MIS mis;
	MES mes;

	for(int i = 0; i < num_vertices(); i++)
	{
		edge_descriptor x = partners[i].first;
		if(x == null_edge) continue;
		assert(mes.find(x) == mes.end());
		mes.insert(PES(x, tostring(i)));
	}

	for(int i = 0; i < num_vertices(); i++)
	{
		edge_descriptor x = partners[i].second;
		if(x == null_edge) continue;
		if(mes.find(x) == mes.end()) mes.insert(PES(x, tostring(i)));
		else mes[x].append("," + tostring(i));
	}
	directed_graph::draw(file, mis, mes, 3.5);
	return 0;
}

