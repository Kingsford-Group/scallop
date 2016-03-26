#include "scallop.h"
#include "subsetsum.h"
#include "nested_graph.h"
#include <cstdio>
#include <iostream>
#include <cfloat>

scallop::scallop(const string &name, splice_graph &gr)
	: assembler(name, gr)
{}

scallop::~scallop()
{}

int scallop::assemble()
{
	smooth_weights();
	init_super_edges();
	reconstruct_splice_graph();
	gr.get_edge_indices(i2e, e2i);
	init_disjoint_sets();
	nt.build(gr);

	round = 0;
	printf("\nprocess %s\n", name.c_str());
	print();

	while(iterate());

	int p = gr.compute_decomp_paths();

	collect_paths();

	printf("%s solution %lu paths, expected = %d paths\n", name.c_str(), paths.size(), p);

	return 0;
}

bool scallop::iterate()
{
	while(true)
	{
		int ei;
		vector<int> sub;
		bool flag0 = identify_equation(ei, sub);
		if(flag0 == true) 
		{
			split_edge(ei, sub);
			nt.build(gr);
			if(sub.size() >= 2)
			{
				printf("decompose %d according to = (#", ei);
				for(int i = 0; i < sub.size(); i++) printf(", %d", sub[i]);
				printf(")\n");

				print();
			}
		}

		bool flag1 = false;
		while(true)
		{
			int ex, ey;
			vector<PI> p;
			bool b1 = identify_linkable_edges(ex, ey, p);
			if(b1 == true)
			{
				printf("linkable edges = (%d, %d), path = (#", ex, ey);
				for(int i = 0; i < p.size(); i++) printf(", (%d,%d)", p[i].first, p[i].second);
				printf(")\n");

				build_adjacent_equal_edges(p);
				connect_adjacent_equal_edges(ex, ey);
				nt.build(gr);

				print();
			}

			bool b2 = decompose_trivial_vertices();
			nt.build(gr);
			if(b2 == true) print();

			if(b1 == true || b2 == true) flag1 = true;
			if(b1 == false && b2 == false) break;
		}

		if(flag0 == false && flag1 == false) break;
	}

	int ex, ey;
	bool b = compute_closest_equal_edges(ex, ey);

	if(b == false) return false;

	printf("shortest equal path (%d, %d)\n", ex, ey);
	connect_equal_edges(ex, ey);
	nt.build(gr);

	print();

	return true;
}

int scallop::init_super_edges()
{
	mev.clear();
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		vector<int> v;
		int s = (*it1)->source();
		v.push_back(s);
		mev.insert(PEV(*it1, v));
	}
	return 0;
}

int scallop::reconstruct_splice_graph()
{
	while(true)
	{
		bool flag = false;
		for(int i = 0; i < gr.num_vertices(); i++)
		{
			bool b = init_trivial_vertex(i);
			if(b == true) flag = true;
		}
		if(flag == false) break;
	}
	return 0;
}

bool scallop::init_trivial_vertex(int x)
{
	int id = gr.in_degree(x);
	int od = gr.out_degree(x);

	if(id <= 0 || od <= 0) return false;
	if(id >= 2 && od >= 2) return false;
	//if(id <= 1 && od <= 1) return false;

	edge_iterator it1, it2;
	edge_iterator ot1, ot2;
	for(tie(it1, it2) = gr.in_edges(x); it1 != it2; it1++)
	{

		for(tie(ot1, ot2) = gr.out_edges(x); ot1 != ot2; ot1++)
		{
			int s = (*it1)->source();
			int t = (*ot1)->target();

			double w1 = gr.get_edge_weight(*it1);
			double a1 = gr.get_edge_stddev(*it1);
			double w2 = gr.get_edge_weight(*ot1);
			double a2 = gr.get_edge_stddev(*ot1);

			double w = w1 < w2 ? w1 : w2;
			double a = w1 < w2 ? a1 : a2;

			edge_descriptor p = gr.add_edge(s, t);
			gr.set_edge_weight(p, w);
			gr.set_edge_stddev(p, a);

			assert(mev.find(*it1) != mev.end());
			assert(mev.find(*ot1) != mev.end());

			vector<int> v1 = mev[*it1];
			vector<int> v2 = mev[*ot1];
			v1.insert(v1.end(), v2.begin(), v2.end());

			if(mev.find(p) != mev.end()) mev[p] = v1;
			else mev.insert(PEV(p, v1));
		}
	}
	gr.clear_vertex(x);
	return true;
}

int scallop::init_disjoint_sets()
{
	ds = disjoint_sets_t(gr.num_edges() * gr.num_vertices());
	for(int i = 0; i < gr.num_edges(); i++)
	{
		ds.make_set(i);
	}
	return 0;
}

vector<int> scallop::compute_representatives()
{
	vector<int> v;
	vector< vector<int> > vv = get_disjoint_sets(ds, i2e.size());
	for(int i = 0; i < vv.size(); i++)
	{
		if(vv[i].size() <= 0) continue;
		int k = -1;
		for(int j = 0; j < vv[i].size(); j++)
		{
			int e = vv[i][j];
			if(i2e[e] == null_edge) continue;
			k = e;
			break;
		}
		if(k == -1) continue;
		v.push_back(k);
	}
	return v;
}

vector< vector<int> > scallop::compute_disjoint_sets()
{
	vector< vector<int> > xx;
	vector< vector<int> > vv = get_disjoint_sets(ds, i2e.size());
	for(int i = 0; i < vv.size(); i++)
	{
		if(vv[i].size() == 0) continue;
		vector<int> v;
		for(int j = 0; j < vv[i].size(); j++)
		{
			int e = vv[i][j];
			if(i2e[e] == null_edge) continue;
			v.push_back(e);
		}
		if(v.size() <= 0) continue;
		xx.push_back(v);
	}
	return xx;
}

set<int> scallop::compute_singletons()
{
	set<int> s;
	vector< vector<int> > vv = compute_disjoint_sets();
	for(int i = 0; i < vv.size(); i++)
	{
		assert(vv[i].size() >= 1);
		if(vv[i].size() >= 2) continue;
		s.insert(vv[i][0]);
	}
	return s;
}

int scallop::connect_equal_edges(int x, int y)
{
	assert(i2e[x] != null_edge);
	assert(i2e[y] != null_edge);

	edge_descriptor xx = i2e[x];
	edge_descriptor yy = i2e[y];

	int xs = (xx)->source();
	int xt = (xx)->target();
	int ys = (yy)->source();
	int yt = (yy)->target();

	double wx = gr.get_edge_weight(xx);
	double wy = gr.get_edge_weight(yy);
	assert(fabs(wx - wy) <= SMIN);

	VE p;
	int l = gr.compute_shortest_path_w(xt, ys, wx, p);
	assert(l >= 0);

	p.insert(p.begin(), xx);
	p.insert(p.end(), yy);

	return connect_path(p, wx);
}

int scallop::connect_path(const VE &p, double wx)
{
	vector<int> v;
	for(int i = 0; i < p.size(); i++)
	{
		assert(p[i] != null_edge);
		assert(e2i.find(p[i]) != e2i.end());
		v.push_back(e2i[p[i]]);
	}
	return connect_path(v, wx);
}

int scallop::connect_adjacent_equal_edges(int x, int y)
{
	if(i2e[x] == null_edge) return -1;
	if(i2e[y] == null_edge) return -1;

	edge_descriptor xx = i2e[x];
	edge_descriptor yy = i2e[y];

	int xs = (xx)->source();
	int xt = (xx)->target();
	int ys = (yy)->source();
	int yt = (yy)->target();

	if(xt != ys && yt != xs) return -1;
	if(yt == xs) return connect_adjacent_equal_edges(y, x);
	
	assert(xt == ys);

	edge_descriptor p = gr.add_edge(xs, yt);

	int n = i2e.size();
	i2e.push_back(p);
	assert(e2i.find(p) == e2i.end());
	e2i.insert(PEI(p, n));

	double wx0 = gr.get_edge_weight(xx);
	double wy0 = gr.get_edge_weight(yy);
	double wx1 = gr.get_edge_stddev(xx);
	double wy1 = gr.get_edge_stddev(yy);

	assert(fabs(wx0 - wy0) <= SMIN);

	gr.set_edge_weight(p, wx0);
	gr.set_edge_stddev(p, wx1);

	vector<int> v = mev[xx];
	v.insert(v.end(), mev[yy].begin(), mev[yy].end());

	if(mev.find(p) != mev.end()) mev[p] = v;
	else mev.insert(PEV(p, v));

	assert(i2e[n] == p);
	assert(e2i.find(p) != e2i.end());
	assert(e2i[p] == n);
	assert(e2i[i2e[n]] == n);

	ds.make_set(n);
	ds.union_set(n, x);
	ds.union_set(n, y);

	e2i.erase(xx);
	e2i.erase(yy);
	i2e[x] = null_edge;
	i2e[y] = null_edge;
	gr.remove_edge(xx);
	gr.remove_edge(yy);

	return n;
}

int scallop::connect_path(const vector<int> &p, double ww)
{
	if(p.size() == 0) return -1;
	if(p.size() == 1) return p[0];

	int pret = i2e[p[0]]->target();
	for(int i = 1; i < p.size(); i++)
	{
		int s = i2e[p[i]]->source();
		int t = i2e[p[i]]->target();
		assert(s == pret);
		pret = t;
	}

	int pree = split_edge(p[0], ww);
	for(int i = 1; i < p.size(); i++)
	{
		int e = split_edge(p[i], pree);
		int ee = connect_adjacent_equal_edges(pree, e);
		pree = ee;
	}
	return pree;
}

int scallop::split_edge(int ei, double w)
{
	assert(i2e[ei] != null_edge);
	edge_descriptor ee = i2e[ei];

	double ww = gr.get_edge_weight(ee);
	double dd = gr.get_edge_stddev(ee);

	if(fabs(ww - w) <= SMIN) return ei;

	assert(ww >= w + SMIN);

	int s = ee->source();
	int t = ee->target();

	edge_descriptor p1 = gr.add_edge(s, t);
	edge_descriptor p2 = gr.add_edge(s, t);

	gr.set_edge_weight(p1, w);
	gr.set_edge_weight(p2, ww - w);
	gr.set_edge_stddev(p1, dd);
	gr.set_edge_stddev(p2, dd);

	if(mev.find(p1) != mev.end()) mev[p1] = mev[ee];
	else mev.insert(PEV(p1, mev[ee]));

	if(mev.find(p2) != mev.end()) mev[p2] = mev[ee];
	else mev.insert(PEV(p2, mev[ee]));

	int n = i2e.size();
	i2e.push_back(p1);
	i2e.push_back(p2);
	e2i.insert(PEI(p1, n));
	e2i.insert(PEI(p2, n + 1));

	gr.remove_edge(ee);
	e2i.erase(ee);
	i2e[ei] = null_edge;

	return n;
}

int scallop::split_edge(int exi, int eyi)
{
	assert(i2e[exi] != null_edge);
	assert(i2e[eyi] != null_edge);
	edge_descriptor ex = i2e[exi];
	edge_descriptor ey = i2e[eyi];

	double wx = gr.get_edge_weight(ex);
	double wy = gr.get_edge_weight(ey);

	int ee = split_edge(exi, wy);

	if(ee == exi)
	{
		ds.union_set(exi, eyi);
		return exi;
	}

	ds.make_set(ee);
	ds.make_set(ee + 1);
	ds.union_set(ee, eyi);

	return ee;
}

vector<int> scallop::split_edge(int ei, const vector<int> &sub)
{
	vector<int> v;
	int x = ei;
	for(int i = 0; i < sub.size(); i++)
	{
		int y = split_edge(x, sub[i]);
		if(i == sub.size() - 1) assert(y == x);
		if(y == x) assert(i == sub.size() - 1);
		v.push_back(y);
		x = y + 1;
	}
	assert(v.size() == sub.size());
	return v;
}

bool scallop::identify_equation(int &ei, vector<int> &sub)
{
	vector<PI> p;
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		int w = (int)(gr.get_edge_weight(i2e[i]));
		p.push_back(PI(w, i));
	}
	sort(p.begin(), p.end());

	sub.clear();
	for(int i = 0; i < p.size(); i++)
	{
		int e = p[i].second;
		assert(i2e[e] != null_edge);
		vector<int> v;
		bool b = identify_edge_equation(e, v);
		if(b == false) continue;
		if(sub.size() == 0 || v.size() < sub.size())
		{
			ei = e;
			sub = v;
		}
	}
	if(sub.size() == 0) return false;
	else return true;
}

bool scallop::identify_edge_equation(int ei, vector<int> &sub)
{
	int s = i2e[ei]->source();
	int t = i2e[ei]->target();
	int w = (int)(gr.get_edge_weight(i2e[ei]));

	set<edge_descriptor> ff;
	set<edge_descriptor> bb;
	gr.bfs(t, ff);
	gr.bfs_reverse(s, bb);

	set<int> sr = compute_singletons();
	//sr.insert(ds.find_set(ei));
	vector<PI> xi;
	for(set<edge_descriptor>::iterator it = ff.begin(); it != ff.end(); it++)
	{
		/*
		int r = ds.find_set(e2i[*it]);
		if(sr.find(r) != sr.end()) continue;
		sr.insert(r);
		*/
		if(sr.find(e2i[*it]) == sr.end()) continue;
		int ww = (int)(gr.get_edge_weight(*it));
		if(ww > w) continue;
		xi.push_back(PI(ww, e2i[*it]));
	}
	for(set<edge_descriptor>::iterator it = bb.begin(); it != bb.end(); it++)
	{
		/*
		int r = ds.find_set(e2i[*it]);
		if(sr.find(r) != sr.end()) continue;
		sr.insert(r);
		*/
		if(sr.find(e2i[*it]) == sr.end()) continue;
		int ww = (int)(gr.get_edge_weight(*it));
		if(ww > w) continue;
		xi.push_back(PI(ww, e2i[*it]));
	}

	sort(xi.begin(), xi.end());

	xi.push_back(PI(w, ei));

	/*
	printf("subsetsum: ");
	for(int i = 0; i < xi.size(); i++)
	{
		printf("%d:%d ", xi[i].first, xi[i].second);
	}
	printf("\n");
	*/

	vector<int> v;
	for(int i = 0; i < xi.size(); i++)
	{
		v.push_back(xi[i].first);
	}

	subsetsum sss(v);
	sss.solve();

	sub.clear();
	for(int j = 0; j < sss.subset.size(); j++)
	{
		int k = sss.subset[j];
		sub.push_back(xi[k].second);
	}

	int err = (int)fabs(sss.opt - w);
	if(err >= 1) return false;

	/*
	printf("closest subset for edge %d:%d has %lu edges, error = %d, subset = (", ei, w, sub.size(), err);
	for(int i = 0; i < sub.size() - 1; i++) printf("%d:%.0lf, ", sub[i], gr.get_edge_weight(i2e[sub[i]]));
	printf("%d:%.0lf)\n", sub[sub.size() - 1], gr.get_edge_weight(i2e[sub[sub.size() - 1]]));
	*/
	
	return true;
}

bool scallop::check_linkable(int ex, int ey, vector<PI> &p)
{
	assert(i2e[ex] != null_edge);
	assert(i2e[ey] != null_edge);

	int xs = i2e[ex]->source();
	int xt = i2e[ex]->target();
	int ys = i2e[ey]->source();
	int yt = i2e[ey]->target();

	vector<PI> xp, yp;
	bool b = nt.link(xs, xt, ys, yt, xp, yp);
	if(b == false) return false;

	p = xp;
	p.insert(p.end(), yp.begin(), yp.end());

	return true;
}

bool scallop::identify_linkable_edges(int &ex, int &ey, vector<PI> &p)
{
	ex = ey = -1;
	p.clear();
	vector< vector<int> > vv = compute_disjoint_sets();
	bool flag = false;
	for(int i = 0; i < vv.size(); i++)
	{
		vector<int> &v = vv[i];
		if(v.size() == 1) continue;
		for(int j = 0; j < v.size(); j++)
		{
			for(int k = j + 1; k < v.size(); k++)
			{
				bool b = check_linkable(v[j], v[k], p);
				//printf(" check %d and %d = %c\n", v[j], v[k], b ? 'T' : 'F');
				if(b == false) continue;
				ex = v[j];
				ey = v[k];
				flag = true;
				break;
			}
			if(flag == true) break;
		}
		if(flag == true) break;
	}
	if(flag == false) return false;
	return true;
}

int scallop::build_adjacent_equal_edges(const vector<PI> &p)
{
	for(int i = 0; i < p.size(); i++)
	{
		int x = p[i].first;
		int y = p[i].second;
		if(y == -1)
		{
			int l = gr.compute_in_partner(x);
			int r = gr.compute_out_partner(x);
			gr.exchange(l, x, r);
		}
		else
		{
			gr.rotate(x, y);
		}
	}
	return 0;
}

bool scallop::decompose_trivial_vertices()
{
	bool flag = false;
	edge_iterator it1, it2;
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		if(gr.degree(i) == 0) continue;
		if(gr.in_degree(i) == 1)
		{
			printf("decompose trivial vertex %d\n", i);

			tie(it1, it2) = gr.in_edges(i);
			int ei = e2i[*it1];
			vector<int> sub;
			for(tie(it1, it2) = gr.out_edges(i); it1 != it2; it1++)
			{
				int e = e2i[*it1];
				sub.push_back(e);
			}

			vector<int> v = split_edge(ei, sub);
			assert(v.size() == sub.size());
			for(int k = 0; k < v.size(); k++)
			{
				assert(i2e[v[k]] != null_edge);
				connect_adjacent_equal_edges(v[k], sub[k]);
			}
			flag = true;
		}
		else if(gr.out_degree(i) == 1)
		{
			printf("decompose trivial vertex %d\n", i);

			tie(it1, it2) = gr.out_edges(i);
			int ei = e2i[*it1];
			vector<int> sub;
			for(tie(it1, it2) = gr.in_edges(i); it1 != it2; it1++)
			{
				int e = e2i[*it1];
				sub.push_back(e);
			}

			vector<int> v = split_edge(ei, sub);
			assert(v.size() == sub.size());
			for(int k = 0; k < v.size(); k++)
			{
				assert(i2e[v[k]] != null_edge);
				connect_adjacent_equal_edges(sub[k], v[k]);
			}
			flag = true;
		}
	}
	return flag;
}

bool scallop::compute_closest_equal_edges(int &ex, int &ey)
{
	ex = ey = -1;
	vector< vector<int> > vv = compute_disjoint_sets();
	int min = INT_MAX;
	for(int i = 0; i < vv.size(); i++)
	{
		vector<int> &v = vv[i];
		if(v.size() == 1) continue;
		for(int j = 0; j < v.size(); j++)
		{
			int xs = i2e[v[j]]->source();
			int xt = i2e[v[j]]->target();
			double wx = gr.get_edge_weight(i2e[v[j]]);
			for(int k = j + 1; k < v.size(); k++)
			{
				int ys = i2e[v[k]]->source();
				int yt = i2e[v[k]]->target();
				double wy = gr.get_edge_weight(i2e[v[k]]);

				assert(fabs(wx - wy) <= SMIN);

				int pxy = gr.compute_shortest_path_w(xt, ys, wx);
				int pyx = gr.compute_shortest_path_w(yt, xs, wx);
				
				assert(pxy == -1 || pyx == -1);
				if(pxy == -1 && pyx == -1) continue;

				if(pxy >= 0 && pxy < min)
				{
					min = pxy;
					ex = v[j];
					ey = v[k];
				}
				else if(pyx >= 0 && pyx < min)
				{
					min = pyx;
					ex = v[k];
					ey = v[j];
				}
			}
		}
	}
	if(min == INT_MAX) return false;
	else return true;
}

int scallop::collect_paths()
{
	paths.clear();
	while(true)
	{
		VE v;
		double w = gr.compute_maximum_path_w(v);
		if(w <= 0.5) break;

		int e = connect_path(v, w);

		path p;
		p.abd = w;
		p.v = recover_path(e);
		paths.push_back(p);

		gr.remove_edge(i2e[e]);
		e2i.erase(i2e[e]);
		i2e[e] = null_edge;
	}
	return 0;
}

vector<int> scallop::recover_path(int e)
{
	assert(i2e[e] != null_edge);
	assert(mev.find(i2e[e]) != mev.end());
	vector<int> v = mev[i2e[e]];
	sort(v.begin(), v.end());

	assert(v[0] == 0);
	assert(v[v.size() - 1] < gr.num_vertices() - 1);
	v.push_back(gr.num_vertices() - 1);
	return v;
}

int scallop::print()
{
	/*
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		vector<int> v = mev[*it1];
		printf("vertices in edge %d = ", e2i[*it1]);
		for(int i = 0; i < v.size(); i++) printf("%d ", v[i]);
		printf("\n");
	}

	vector< vector<int> > vv = compute_disjoint_sets();
	for(int i = 0; i < vv.size(); i++)
	{
		vector<int> v = vv[i];
		assert(v.size() >= 1);

		if(v.size() == 1) continue;
		int w = (int)(gr.get_edge_weight(i2e[v[0]]));

		printf("edge set %d, weight = %d, #edges = %lu, set = (%d", i, w, v.size(), v[0]);
		for(int j = 1; j < v.size(); j++) printf(", %d", v[j]);
		printf(")\n");
	}
	*/
	int n = 0;
	for(int i = 0; i < gr.num_vertices(); i++) 
	{
		if(gr.degree(i) >= 1) n++;
	}
	printf("statistics: %lu edges, %d vertices\n", gr.num_edges(), n);
	printf("finish round %d\n\n", round);

	char buf[1024];
	sprintf(buf, "%s.gr.%d.tex", name.c_str(), round);
	draw_splice_graph(buf);

	sprintf(buf, "%s.nt.%d.tex", name.c_str(), round);
	nt.draw(buf);

	round++;

	return 0;
}

int scallop::draw_splice_graph(const string &file) 
{
	MIS mis;
	char buf[10240];
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		double w = gr.get_vertex_weight(i);
		sprintf(buf, "%d:%.0lf", i, w);
		mis.insert(PIS(i, buf));
	}

	MES mes;
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		double w = gr.get_edge_weight(i2e[i]);
		sprintf(buf, "%d:%.0lf", i, w);
		//sprintf(buf, "%d", i);
		mes.insert(PES(i2e[i], buf));
	}
	gr.draw(file, mis, mes, 4.5);
	return 0;
}
