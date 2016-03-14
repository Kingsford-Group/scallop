#include "scallop3.h"
#include "subsetsum.h"
#include <cstdio>
#include <iostream>

using namespace subsetsum;

scallop3::scallop3(const string &name, splice_graph &gr)
	: assembler(name, gr)
{}

scallop3::~scallop3()
{}

int scallop3::assemble()
{
	smooth_weights();
	init_super_edges();
	reconstruct_splice_graph();
	get_edge_indices(gr, i2e, e2i);
	init_disjoint_sets();
	build_null_space();
	iterate();
	return 0;
}

int scallop3::iterate()
{
	while(true)
	{
		int ei;
		vector<int> sub;
		int error = identify_equation(ei, sub);
		if(error >= 1) break;
		bool b = verify_equation(ei, sub);
		if(b == false) break;
		split_edge(ei, sub);
	}

	this->draw_splice_graph(name + ".gr.tex");

	while(true)
	{
		build_nested_graph();
		int ex, ey;
		vector<int> p;
		bool b = identify_linkable_edges(ex, ey, p);
		if(b == false) break;
		assert(p.size() >= 1);
		printf("linkable edges = (%d, %d), path = (%d", ex, ey, p[0]);
		for(int i = 1; i < p.size(); i++) printf(", %d", p[i]);
		printf(")\n");
		assert(ex >= 0 && ey >= 0);
		//connect_adjacent_edges(ex, ey);
		break;
	}

	this->draw_nested_graph(name + ".nt.tex");

	return 0;
}

int scallop3::init_super_edges()
{
	mev.clear();
	edge_iterator it1, it2;
	vector<int> v;
	for(tie(it1, it2) = edges(gr); it1 != it2; it1++)
	{
		mev.insert(PEV(*it1, v));
	}
	return 0;
}

int scallop3::reconstruct_splice_graph()
{
	while(true)
	{
		bool flag = false;
		for(int i = 0; i < num_vertices(gr); i++)
		{
			bool b = decompose_trivial_vertex(i);
			if(b == true) flag = true;
		}
		if(flag == false) break;
	}
	return 0;
}

bool scallop3::decompose_trivial_vertex(int x)
{
	int id = in_degree(x, gr);
	int od = out_degree(x, gr);

	if(id <= 0 || od <= 0) return false;
	if(id >= 2 && od >= 2) return false;
	//if(id <= 1 && od <= 1) return false;

	in_edge_iterator it1, it2;
	out_edge_iterator ot1, ot2;
	for(tie(it1, it2) = in_edges(x, gr); it1 != it2; it1++)
	{

		for(tie(ot1, ot2) = out_edges(x, gr); ot1 != ot2; ot1++)
		{
			int s = source(*it1, gr);
			int t = target(*ot1, gr);

			double w1 = get_edge_weight(*it1, gr);
			double a1 = get_edge_stddev(*it1, gr);
			double w2 = get_edge_weight(*ot1, gr);
			double a2 = get_edge_stddev(*ot1, gr);

			double w = w1 < w2 ? w1 : w2;
			double a = w1 < w2 ? a1 : a2;

			PEB p = add_edge(s, t, gr);
			set_edge_weight(p.first, w, gr);
			set_edge_stddev(p.first, a, gr);

			assert(mev.find(*it1) != mev.end());
			assert(mev.find(*ot1) != mev.end());

			vector<int> v = mev[*it1];
			v.push_back(x);
			v.insert(v.end(), mev[*ot1].begin(), mev[*ot1].end());

			mev.insert(PEV(p.first, v));
		}
	}
	clear_vertex(x, gr);
	return true;
}

int scallop3::init_disjoint_sets()
{
	ds = disjoint_sets_t(num_edges(gr) * num_vertices(gr));
	for(int i = 0; i < num_edges(gr); i++)
	{
		ds.make_set(i);
	}
	return 0;
}

int scallop3::build_null_space()
{
	ns.clear();
	for(int i = 1; i < num_vertices(gr) - 1; i++)
	{
		if(degree(i, gr) == 0) continue;
		vector<double> v;
		v.resize(num_edges(gr), 0);
		in_edge_iterator i1, i2;
		for(tie(i1, i2) = in_edges(i, gr); i1 != i2; i1++)
		{
			int e = e2i[*i1];
			v[e] = 1;
		}
		out_edge_iterator o1, o2;
		for(tie(o1, o2) = out_edges(i, gr); o1 != o2; o1++)
		{
			int e = e2i[*o1];
			v[e] = -1;
		}
		ns.push_back(v);
	}

	algebra ab(ns);
	int rank = ab.partial_eliminate();
	assert(rank == ns.size());

	return 0;
}

vector<int> scallop3::compute_representatives()
{
	vector<int> v;
	vector< vector<int> > vv = get_disjoint_sets(ds, i2e.size());
	for(int i = 0; i < vv.size(); i++)
	{
		if(vv[i].size() == 0) continue;
		int k = -1;
		for(int j = 0; j < vv[i].size(); j++)
		{
			int e = vv[i][j];
			if(i2e[e] == null_edge) continue;
			k = e;
			break;
		}
		assert(k != -1);
		v.push_back(k);
	}
	return v;
}

vector< vector<int> > scallop3::compute_disjoint_sets()
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
		assert(v.size() >= 1);
		xx.push_back(v);
	}
	return xx;
}

bool scallop3::connect_adjacent_edges(int x, int y)
{
	if(i2e[x] == null_edge) return false;
	if(i2e[y] == null_edge) return false;

	edge_descriptor xx = i2e[x];
	edge_descriptor yy = i2e[y];

	int xs = source(xx, gr);
	int xt = target(xx, gr);
	int ys = source(yy, gr);
	int yt = target(yy, gr);

	if(xt != ys && yt != xs) return false;
	if(yt == xs) return connect_adjacent_edges(y, x);
	
	assert(xt == ys);

	PEB p = add_edge(xs, yt, gr);
	assert(p.second == true);

	int n = i2e.size();
	i2e.push_back(p.first);
	assert(e2i.find(p.first) == e2i.end());
	e2i.insert(PEI(p.first, n));

	double wx0 = get_edge_weight(xx, gr);
	double wy0 = get_edge_weight(yy, gr);
	double wx1 = get_edge_stddev(xx, gr);
	double wy1 = get_edge_stddev(yy, gr);

	//printf("wx0 = %lf, wy0 = %lf, SMIN = %lf\n", wx0, wy0, SMIN);

	assert(fabs(wx0 - wy0) <= SMIN);

	set_edge_weight(p.first, wx0, gr);
	set_edge_stddev(p.first, wx1, gr);

	vector<int> v = mev[xx];
	v.push_back(xt);
	v.insert(v.end(), mev[yy].begin(), mev[yy].end());

	mev.insert(PEV(p.first, v));

	assert(i2e[n] == p.first);
	assert(e2i.find(p.first) != e2i.end());
	assert(e2i[p.first] == n);
	assert(e2i[i2e[n]] == n);

	ds.make_set(n);
	ds.union_set(n, x);
	ds.union_set(n, y);

	e2i.erase(xx);
	e2i.erase(yy);
	i2e[x] = null_edge;
	i2e[y] = null_edge;
	remove_edge(xx, gr);
	remove_edge(yy, gr);

	return 0;
}

int scallop3::split_edge(int ei, const vector<int> &sub)
{
	assert(i2e[ei] != null_edge);
	assert(sub.size() >= 1);
	for(int i = 0; i < sub.size(); i++) assert(i2e[sub[i]] != null_edge);

	double w = get_edge_weight(i2e[ei], gr);
	for(int i = 0; i < sub.size(); i++)
	{
		w -= get_edge_weight(i2e[sub[i]], gr);
	}
	assert(fabs(w) <= SMIN);

	edge_descriptor ex = i2e[ei];
	edge_descriptor ey = i2e[sub[0]];
	double w0 = get_edge_weight(ey, gr);
	double w1 = get_edge_stddev(ey, gr);
	set_edge_weight(ex, w0, gr);
	set_edge_stddev(ex, w1, gr);
	ds.union_set(ei, sub[0]);

	int s = source(ex, gr);
	int t = target(ex, gr);
	for(int i = 1; i < sub.size(); i++)
	{
		PEB p = add_edge(s, t, gr);
		assert(p.second == true);

		int n = i2e.size();
		i2e.push_back(p.first);
		assert(e2i.find(p.first) == e2i.end());
		e2i.insert(PEI(p.first, n));

		ey = i2e[sub[i]];
		w0 = get_edge_weight(ey, gr);
		w1 = get_edge_stddev(ey, gr);

		set_edge_weight(p.first, w0, gr);
		set_edge_stddev(p.first, w1, gr);

		mev.insert(PEV(p.first, mev[ex]));

		ds.make_set(n);
		ds.union_set(n, sub[i]);
	}

	// TODO, update null space
	return 0;
}

int scallop3::identify_equation(int &ei, vector<int> &sub)
{
	// TODO DEBUG
	//if(num_edges(gr) >= 11) return 0;

	vector<int> r = compute_representatives();

	vector<int> x;
	for(int i = 0; i < r.size(); i++)
	{
		double w = get_edge_weight(i2e[r[i]], gr);
		x.push_back((int)(w));
	}

	vector<int> xx;
	vector<int> xf;
	vector<int> xb;
	enumerate_subsets(x, xx, xf, xb);

	vector<PI> xxp;
	for(int i = 0; i < xx.size(); i++) xxp.push_back(PI(xx[i], i));

	sort(xxp.begin(), xxp.end());

	// sort all edges? TODO
	vector<PI> xp;
	for(int i = 0; i < x.size(); i++) xp.push_back(PI(x[i], i));

	sort(xp.begin(), xp.end());

	int ri = -1;
	int xxpi = -1;
	int minw = INT_MAX;
	for(int i = 0; i < xp.size(); i++)
	{
		int k = compute_closest_subset(xp[i].second, xp[i].first, xxp);
		double ww = (int)fabs(xp[i].first - xxp[k].first);
		if(ww < minw)
		{
			minw = ww;
			xxpi = k;
			ri = xp[i].second;
		}
	}

	assert(ri >= 0 && ri < r.size());

	ei = r[ri];

	vector<int> rsub;
	recover_subset(rsub, xxp[xxpi].second, xf, xb);

	for(int i = 0; i < rsub.size(); i++) sub.push_back(r[rsub[i]]);

	printf("%s closest subset for edge %d:%d has %lu edges, error = %d, subset = (", name.c_str(), ei, x[ri], sub.size(), minw);
	for(int i = 0; i < sub.size() - 1; i++) printf("%d:%d, ", sub[i], x[rsub[i]]);
	printf("%d:%d), total %lu combinations\n", sub[sub.size() - 1], x[rsub[sub.size() - 1]], xx.size());
	
	return minw;
}

bool scallop3::verify_equation(int ei, const vector<int> &sub)
{
	assert(i2e[ei] != null_edge);
	for(int i = 0; i < sub.size(); i++)
	{
		assert(i2e[sub[i]] != null_edge);
		bool b1 = gr.check_directed_path(i2e[ei], i2e[sub[i]]);
		bool b2 = gr.check_directed_path(i2e[sub[i]], i2e[ei]);
		if(b1 == false && b2 == false) return false;
	}
	return true;
}

int scallop3::compute_closest_subset(int xi, int w, const vector<PI> &xxp)
{
	if(xxp.size() == 0) return -1;

	int s = 0; 
	int t = xxp.size() - 1;
	int m = -1;
	while(s < t)
	{
		m = (s + t) / 2;
		if(w == xxp[m].first) break;
		else if(w > xxp[m].first) s = m;
		else t = m;
	}

	int si = -1;
	for(si = m; si >= 0; si--)
	{
		if(xxp[si].second != xi) break;
	}

	int ti = -1;
	for(ti = m + 1; ti < xxp.size(); ti++)
	{
		if(xxp[ti].second != xi) break;
	}

	assert(si != -1 || ti != -1);
	
	if(si == -1) return ti;
	if(ti == -1) return si;

	int sw = (int)(fabs(xxp[si].first - w));
	int tw = (int)(fabs(xxp[ti].first - w));

	if(sw <= tw) return si;
	else return ti;
}

int scallop3::compute_closest_equal_edges(int &ex, int &ey)
{
	int dist = INT_MAX;
	ex = ey = -1;
	vector< vector<int> > vv = compute_disjoint_sets();
	for(int i = 0; i < vv.size(); i++)
	{
		vector<int> ss;
		vector<int> tt;
		for(int j = 0; j < vv[i].size(); j++)
		{
			int e = vv[i][j];
			int s = source(i2e[e], gr);
			int t = target(i2e[e], gr);
			ss.push_back(s);
			tt.push_back(t);
		}

		for(int j = 0; j < tt.size(); j++)
		{
			for(int k = 0; k < ss.size(); k++)
			{
				if(tt[j] != ss[k]) continue;
				dist = 0;
				ex = vv[i][j];
				ey = vv[i][k];
				return dist;
			}
		}

		for(int j = 0; j < tt.size(); j++)
		{
			vector<int> v;
			bfs_splice_graph(gr, tt[j], v);
			for(int k = 0; k < ss.size(); k++)
			{
				if(v[ss[k]] <= 0) continue;
				if(v[ss[k]] < dist)
				{
					dist = v[ss[k]];
					ex = vv[i][j];
					ey = vv[i][k];
				}
			}
		}
	}
	return dist;
}

int scallop3::build_nested_graph()
{
	nt.clear();
	i2n.clear();
	for(int i = 0; i < gr.num_vertices(); i++) nt.add_vertex();

	i2n.assign(i2e.size(), null_edge);
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		int s = i2e[i]->source();
		int t = i2e[i]->target();
		edge_descriptor e = nt.add_edge(s, t);
		i2n[i] = e;
	}
	return 0;
}

bool scallop3::identify_linkable_edges(int &ex, int &ey, vector<int> &p)
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
				edge_descriptor ej = i2n[v[j]];
				edge_descriptor ek = i2n[v[k]];
				if(ej == null_edge) continue;
				if(ek == null_edge) continue;
				bool b1 = nt.check_directed_path(ej, ek);
				bool b2 = nt.check_directed_path(ek, ej);
				if(b1 == false && b2 == false) continue;
				if(b1 == true)
				{
					ex = v[j];
					ey = v[k];
					nt.compute_shortest_path(ej, ek, p);
				}
				else
				{
					ex = v[k];
					ey = v[j];
					nt.compute_shortest_path(ek, ej, p);
				}
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

int scallop3::print()
{
	// print null space
	/*
	if(ns.size() == 0) return 0;
	printf("null space:\n");
	algebra::print_matrix(ns);
	*/

	// print edge disjoint sets
	vector< vector<int> > vv = compute_disjoint_sets();
	for(int i = 0; i < vv.size(); i++)
	{
		vector<int> v = vv[i];
		assert(v.size() >= 1);

		int w = (int)(get_edge_weight(i2e[v[0]], gr));

		printf("edge set %d, weight = %d, #edges = %lu, set = (%d", i, w, v.size(), v[0]);
		for(int j = 1; j < v.size(); j++) printf(", %d", v[j]);
		printf(")\n");
	}
	return 0;
}

int scallop3::draw_splice_graph(const string &file) 
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
		mes.insert(PES(i2e[i], buf));
	}

	gr.draw(file, mis, mes, 5.0);
	return 0;
}

int scallop3::draw_nested_graph(const string &file) 
{
	MIS mis;
	MES mes;
	char buf[1024];
	for(int i = 0; i < i2n.size(); i++)
	{
		if(i2n[i] == null_edge) continue;
		sprintf(buf, "%d", i);
		mes.insert(PES(i2n[i], buf));
	}
	nt.draw(file, mis, mes, 4.0);
	return 0;
}
