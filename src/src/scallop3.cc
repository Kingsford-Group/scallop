#include "scallop3.h"
#include "subsetsum.h"
#include <cstdio>

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
	draw_splice_graph(gr, name + ".tex", 5.0);
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
		process_equation(ei, sub);
		print();
	}
	return 0;
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
	vector< vector<int> > vv = get_disjoint_sets(ds, i2e.size());
	for(int i = 0; i < vv.size(); i++)
	{
		if(vv[i].size() <= 1) continue;

		vector<int> v;
		for(int j = 0; j < vv[i].size(); j++)
		{
			if(e2i[i2e[vv[i][j]]] == -1) continue;
			v.push_back(vv[i][j]);
		}
		assert(v.size() >= 1);
		
		int w = (int)(get(get(edge_weight, gr), i2e[v[0]]));

		printf("edge set %d, weight = %d, edges = (%d", i, w, v[0]);
		for(int j = 1; j < v.size(); j++) printf(", %d", v[j]);
		printf(")\n");
	}
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

			double w1 = get(get(edge_weight, gr), *it1);
			double a1 = get(get(edge_stddev, gr), *it1);
			double w2 = get(get(edge_weight, gr), *ot1);
			double a2 = get(get(edge_stddev, gr), *ot1);

			double w = w1 < w2 ? w1 : w2;
			double a = w1 < w2 ? a1 : a2;

			PEB p = add_edge(s, t, gr);
			put(get(edge_weight, gr), p.first, w);
			put(get(edge_stddev, gr), p.first, a);

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
			if(e2i[i2e[e]] == -1) continue;
			k = e;
			break;
		}
		assert(k != -1);
		v.push_back(k);
	}
	return v;
}

bool scallop3::connect_edges(int x, int y)
{
	return false;

	edge_descriptor xx = i2e[x];
	edge_descriptor yy = i2e[y];

	if(e2i[xx] == -1) return false;
	if(e2i[yy] == -1) return false;

	int xs = source(xx, gr);
	int xt = target(xx, gr);
	int ys = source(yy, gr);
	int yt = target(yy, gr);

	if(xt != ys && yt != xs) return false;
	if(yt == xs) return connect_edges(y, x);
	
	assert(xt == ys);

	PEB p = add_edge(xs, yt, gr);
	assert(p.second == true);

	int n = i2e.size();
	e2i.insert(PEI(p.first, n));
	i2e.push_back(p.first);

	double wx0 = get(get(edge_weight, gr), xx);
	double wy0 = get(get(edge_weight, gr), yy);
	double wx1 = get(get(edge_stddev, gr), xx);
	double wy1 = get(get(edge_stddev, gr), yy);

	assert(fabs(wx0 - wy0) <= SMIN);

	put(get(edge_weight, gr), p.first, wx0);
	put(get(edge_stddev, gr), p.first, wx1);

	vector<int> v = mev[xx];
	v.push_back(xt);
	v.insert(v.end(), mev[yy].begin(), mev[yy].end());

	mev.insert(PEV(p.first, v));

	ds.make_set(n);
	ds.union_set(n, e2i[xx]);

	e2i[xx] = -1;
	e2i[yy] = -1;
	remove_edge(xx, gr);
	remove_edge(yy, gr);

	return 0;
}

int scallop3::process_equation(int ei, const vector<int> &sub)
{
	assert(sub.size() >= 1);
	assert(e2i[i2e[ei]] != -1);
	for(int i = 0; i < sub.size(); i++) assert(e2i[i2e[sub[i]]] != -1);

	edge_descriptor ex = i2e[ei];
	edge_descriptor ey = i2e[sub[0]];
	double w0 = get(get(edge_weight, gr), ey);
	double w1 = get(get(edge_stddev, gr), ey);
	put(get(edge_weight, gr), ex, w0);
	put(get(edge_stddev, gr), ex, w1);
	ds.union_set(ei, sub[0]);
	connect_edges(ei, sub[0]);

	int s = source(ex, gr);
	int t = target(ex, gr);
	for(int i = 1; i < sub.size(); i++)
	{
		PEB p = add_edge(s, t, gr);
		assert(p.second == true);

		int n = i2e.size();
		e2i.insert(PEI(p.first, n));
		i2e.push_back(p.first);

		ey = i2e[sub[i]];
		w0 = get(get(edge_weight, gr), ey);
		w1 = get(get(edge_stddev, gr), ey);

		put(get(edge_weight, gr), p.first, w0);
		put(get(edge_stddev, gr), p.first, w1);

		mev.insert(PEV(p.first, mev[ex]));

		ds.make_set(n);
		ds.union_set(n, sub[i]);

		connect_edges(n, sub[i]);
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
		double w = get(get(edge_weight, gr), i2e[r[i]]);
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
		int k = locate_closest_subset(xp[i].second, xp[i].first, xxp);
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

int scallop3::locate_closest_subset(int xi, int w, const vector<PI> &xxp)
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
