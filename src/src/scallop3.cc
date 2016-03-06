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
	build_null_space();
	iterate();
	return 0;
}

int scallop3::iterate()
{
	int ei;
	vector<int> sub;
	identify_equation(ei, sub);
	return 0;
}

int scallop3::print()
{
	if(ns.size() == 0) return 0;
	printf("null space:\n");
	algebra::print_matrix(ns);
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

int scallop3::identify_equation(int &ei, vector<int> &sub)
{
	if(num_edges(gr) >= 11) return 0;

	vector<int> x;
	for(int i = 0; i < i2e.size(); i++)
	{
		double w = get(get(edge_weight, gr), i2e[i]);
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

	ei = -1;
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
			ei = xp[i].second;
		}
	}

	assert(ei != -1);

	recover_subset(sub, xxp[xxpi].second, xf, xb);

	printf("%s closest subset for edge %d:%d has %lu edges, error = %d, subset = (", name.c_str(), ei, x[ei], sub.size(), minw);
	for(int i = 0; i < sub.size() - 1; i++) printf("%d:%d, ", sub[i], x[sub[i]]);
	printf("%d:%d)\n", sub[sub.size() - 1], x[sub[sub.size() - 1]]);
	
	return 0;
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
