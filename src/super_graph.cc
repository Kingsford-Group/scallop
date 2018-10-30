/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "super_graph.h"
#include "config.h"
#include <algorithm>
#include <cfloat>

super_graph::super_graph(const splice_graph &gr, const hyper_set &hs)
	:root(gr), hyper(hs)
{}

super_graph::~super_graph()
{}

int super_graph::build()
{
	subs.clear();
	hss.clear();

	build_undirected_graph();
	split_splice_graph();

	return 0;
}

int super_graph::build_undirected_graph()
{
	ug.clear();
	for(int i = 0; i < root.num_vertices(); i++)
	{
		ug.add_vertex();
	}

	edge_iterator it1, it2;
	PEEI pei;
	for(pei = root.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int s = e->source();
		int t = e->target();
		if(s == 0) continue;
		if(t == root.num_vertices() - 1) continue;

		ug.add_edge(s, t);
	}

	return 0;
}

int super_graph::split_splice_graph()
{
	vector< set<int> > vv = ug.compute_connected_components();
	a2b.clear();
	b2a.clear();
	int index = 0;
	for(int k = 0; k < vv.size(); k++)
	{
		set<int> &s = vv[k];
		if(s.size() == 1 && *(s.begin()) == 0) continue;
		if(s.size() == 1 && *(s.begin()) == root.num_vertices() - 1) continue;
		splice_graph gr;
		hyper_set hs;
		split_single_splice_graph(gr, hs, s, index);
		gr.chrm = root.chrm;
		gr.strand = root.strand;
		subs.push_back(gr);
		hss.push_back(hs);
		index++;
	}
	return 0;
}

int super_graph::split_single_splice_graph(splice_graph &gr, hyper_set &hs, const set<int> &ss, int index)
{
	//printf("build single splice graph with index = %d\n", index);
	gr.clear();
	vector<int> vv(ss.begin(), ss.end());
	sort(vv.begin(), vv.end());
	assert(vv.size() >= 1);

	int32_t lpos = root.get_vertex_info(vv[0]).lpos;
	int32_t rpos = root.get_vertex_info(vv[vv.size() - 1]).rpos;

	// vertices
	gr.add_vertex();
	vertex_info vi0;
	vi0.lpos = lpos;
	vi0.rpos = lpos;
	gr.set_vertex_weight(0, 0);
	gr.set_vertex_info(0, vi0);

	for(int i = 0; i < vv.size(); i++)
	{
		int k = vv[i];
		gr.add_vertex();
		gr.set_vertex_weight(i + 1, root.get_vertex_weight(k));
		gr.set_vertex_info(i + 1, root.get_vertex_info(k));
		PI p(index, i + 1);
		a2b.insert(pair<int, PI>(k, p));
		b2a.insert(pair<PI, int>(p, k));
	}

	gr.add_vertex();
	vertex_info vin;
	vin.lpos = rpos;
	vin.rpos = rpos;
	gr.set_vertex_weight(vv.size() + 1, 0);
	gr.set_vertex_info(vv.size() + 1, vin);

	// edges
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = root.out_edges(0), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int t = (*it1)->target();
		if(ss.find(t) == ss.end()) continue;
		assert(a2b.find(t) != a2b.end());
		assert(a2b[t].first == index);
		int y = a2b[t].second;

		edge_descriptor e = gr.add_edge(0, y);
		gr.set_edge_weight(e, root.get_edge_weight(*it1));
		gr.set_edge_info(e, root.get_edge_info(*it1));
	}

	int n = root.num_vertices() - 1;
	for(int i = 0; i < vv.size(); i++)
	{
		int s = vv[i];
		assert(s != 0 && s != n);
		assert(a2b.find(s) != a2b.end());
		assert(a2b[s].first == index);
		int x = a2b[s].second;

		for(pei = root.out_edges(s), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			int t = (*it1)->target();
			assert(t == n || ss.find(t) != ss.end());
			assert(t == n || a2b.find(t) != a2b.end());
			assert(t == n || a2b[t].first == index);
			int y = ((t == n) ? gr.num_vertices() - 1 : a2b[t].second);

			edge_descriptor e = gr.add_edge(x, y);
			gr.set_edge_weight(e, root.get_edge_weight(*it1));
			gr.set_edge_info(e, root.get_edge_info(*it1));
		}
	}

	// hyper-set
	hs.clear();
	for(MVII::iterator it = hyper.nodes.begin(); it != hyper.nodes.end(); it++)
	{
		vector<int> v = it->first;
		int c = it->second;

		bool b = true;
		vector<int> vv;
		for(int k = 0; k < v.size(); k++)
		{
			if(ss.find(v[k]) == ss.end()) b = false;
			if(b == false) break;
			assert(a2b.find(v[k]) != a2b.end());
			assert(a2b[v[k]].first == index);
			int x = a2b[v[k]].second;
			vv.push_back(x);
		}

		if(b == false) continue;

		for(int i = 0; i < vv.size(); i++) vv[i]--;
		hs.add_node_list(vv, c);
	}

	return 0;
}

int super_graph::get_root_vertex(int s, int x) const
{
	assert(s >= 0 && s < subs.size());
	assert(x >= 0 && x < subs[s].num_vertices());
	if(x == 0) return 0;
	if(x == subs[s].num_vertices() - 1) return root.num_vertices() - 1;
	PI p(s, x);
	map<PI, int>::const_iterator it = b2a.find(p);
	assert(it != b2a.end());
	int v = it->second;
	assert(v > 0 && v < root.num_vertices() - 1);
	return v;
}

vector<int> super_graph::get_root_vertices(int s, const vector<int> &x) const
{
	vector<int> vv;
	for(int i = 0; i < x.size(); i++)
	{
		int v = get_root_vertex(s, x[i]);
		vv.push_back(v);
	}
	return vv;
}

bool super_graph::cut_splice_graph()
{
	bool flag = false;
	for(int i = 0; i < subs.size(); i++)
	{
		bool b = cut_single_splice_graph(subs[i], i);
		if(b == true) flag = true;
	}
	return flag;
}

bool super_graph::cut_single_splice_graph(splice_graph &gr, int index)
{
	vector<SE> cuts;
	cuts.resize(gr.num_vertices() - 1);
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int s = e->source();
		int t = e->target();
		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;
		double w = gr.get_edge_weight(e);
		for(int k = s; k < t; k++) cuts[k].insert(e);
	}

	vector<int> ss, tt;
	for(pei = gr.out_edges(0), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int t = (*it1)->target();
		ss.push_back(t);
	}
	for(pei = gr.in_edges(gr.num_vertices() - 1), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		tt.push_back(s);
	}

	double max_sum = 3;
	int min_size = 5;

	double ksum = max_sum + 1.0, kave = 0, kmin = 0;
	int ks = -1, kt = -1;
	VE ke;
	for(int i = 0; i < ss.size(); i++)
	{
		int s = ss[i];
		for(int j = 0; j < tt.size(); j++)
		{
			int t = tt[j];
			if(s >= t) continue;
			if(t - s + 1 < min_size) continue;

			SE &sse = cuts[s - 1];
			SE &tte = cuts[t];
			int cnt1 = 0, cnt2 = 0;
			double sum1 = 0, sum2 = 0;
			VE v;
			for(SE::iterator it = sse.begin(); it != sse.end(); it++)
			{
				edge_descriptor e = (*it);
				if(tte.find(e) != tte.end()) continue;
				v.push_back(e);
				cnt1++;
				sum1 += gr.get_edge_weight(e);	
			}
			for(SE::iterator it = tte.begin(); it != tte.end(); it++)
			{
				edge_descriptor e = (*it);
				if(sse.find(e) != sse.end()) continue;
				v.push_back(e);
				cnt2++;
				sum2 += gr.get_edge_weight(e);
			}

			//if(cnt1 >= 1 && s < min_size) continue;
			//if(cnt2 >= 1 && gr.num_vertices() - 1 - t < min_size) continue;
			if(sum1 + sum2 > max_sum) continue;
			if(s <= 1 && t >= gr.num_vertices() - 2) continue;

			double ave = 0, min = DBL_MAX;
			for(int k = s; k < t; k++)
			{
				SE &se = cuts[k];
				double w = 0;
				for(SE::iterator it = se.begin(); it != se.end(); it++)
				{
					edge_descriptor e = (*it);
					if(sse.find(e) != sse.end()) continue;
					if(tte.find(e) != tte.end()) continue;
					w += gr.get_edge_weight(e);
				}
				ave += w;
				if(w < min) min = w;
			}
			ave /= (t - s);
			if(3.0 * (sum1 + sum2) >= ave) continue;

			/*
			int32_t p1 = gr.get_vertex_info(s).lpos;
			int32_t p2 = gr.get_vertex_info(t).rpos;
			printf(" cuts [%d, %d] out of %lu, cnt = (%d, %d), sum = (%.0lf, %.0lf), pos = %d-%d\n", s, t, gr.num_vertices(), cnt1, cnt2, sum1, sum2, p1, p2);
			*/

			if(sum1 + sum2 < ksum)
			{
				ks = s;
				kt = t;
				ke = v;
				ksum = sum1 + sum2;
				kave = ave;
				kmin = min;
			}
		}
	}

	if(ks == -1 || kt == -1) return false;

	printf("cut subgraph %d, vertices = [%d, %d] / %lu, #edges = %.0lf, ave = %.2lf, min = %.2lf\n", index, ks, kt, gr.num_vertices(), ksum, kave, kmin);

	for(int i = 0; i < ke.size(); i++)
	{
		edge_descriptor e = ke[i];
		int s = e->source();
		int t = e->target();
		int ss = b2a[PI(index, s)];
		int tt = b2a[PI(index, t)];
		PEB p = root.edge(ss, tt);
		assert(p.second == true);
		root.remove_edge(p.first);
	}

	return true;
}

int super_graph::build_maximum_path_graph(splice_graph &gr, undirected_graph &mg)
{
	mg.clear();
	int n = gr.num_vertices() - 1;
	if(gr.out_degree(0) <= 1) return 0;
	if(gr.in_degree(n) <= 1) return 0;

	edge_iterator it1, it2;
	PEEI pei;
	vector<int> u2v;
	MI v2u;
	for(pei = gr.out_edges(0), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int s = e->target();
		assert(v2u.find(s) == v2u.end());
		v2u.insert(PI(s, v2u.size()));
		u2v.push_back(s);
		mg.add_vertex();
	}
	for(pei = gr.in_edges(n), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int t = e->source();
		assert(v2u.find(t) == v2u.end());
		v2u.insert(PI(t, v2u.size()));
		u2v.push_back(t);
		mg.add_vertex();
	}

	for(int i = 0; i < u2v.size(); i++)
	{
		int x = u2v[i];
		int y = -1;
		if(i < gr.out_degree(0)) 
		{
			compute_maximum_path1(gr, x, y);
			int j = v2u[y];
			assert(j >= gr.out_degree(0));
			assert(j < gr.out_degree(0) + gr.in_degree(n));
			mg.add_edge(i, j);
		}
		else if(i >= gr.out_degree(0) && i < gr.out_degree(0) + gr.in_degree(n))
		{
			compute_maximum_path2(gr, x, y);
			int j = v2u[y];
			assert(j >= 0);
			assert(j < gr.out_degree(0));
			mg.add_edge(j, i);
		}
		else assert(false);
	}

	return 0;
}

double super_graph::compute_maximum_path1(splice_graph &gr, int s, int &t)
{
	MED med;
	int n = gr.num_vertices() - 1;
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.in_edges(n), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		double w = gr.get_edge_weight(e);
		med.insert(PED(e, w));
		gr.set_edge_weight(e, DBL_MAX / 10.0);
	}

	VE ve;
	double ww = gr.compute_maximum_st_path_w(ve, s, n);
	assert(ve.size() >= 1);
	edge_descriptor e = ve[ve.size() - 1];
	assert(e->target() == n);
	t = e->source();

	for(MED::iterator it = med.begin(); it != med.end(); it++)
	{
		edge_descriptor e = it->first;
		double w = it->second;
		gr.set_edge_weight(e, w);
	}

	return ww;
}

double super_graph::compute_maximum_path2(splice_graph &gr, int t, int &s)
{
	MED med;
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.out_edges(0), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		double w = gr.get_edge_weight(e);
		med.insert(PED(e, w));
		gr.set_edge_weight(e, DBL_MAX / 10.0);
	}

	VE ve;
	double ww = gr.compute_maximum_st_path_w(ve, 0, t);
	assert(ve.size() >= 1);
	edge_descriptor e = ve[0];
	assert(e->source() == 0);
	s = e->target();

	for(MED::iterator it = med.begin(); it != med.end(); it++)
	{
		edge_descriptor e = it->first;
		double w = it->second;
		gr.set_edge_weight(e, w);
	}

	return ww;
}

int super_graph::print()
{
	int d0 = root.out_degree(0);
	int dn = root.in_degree(root.num_vertices() - 1);
	printf("super graph, %lu vertices, starting = %d, ending = %d, %lu subgraphs\n",
			root.num_vertices(), d0, dn, subs.size());

	for(int i = 0; i < subs.size(); i++)
	{
		splice_graph &gr = subs[i];
		d0 = gr.out_degree(0);
		dn = gr.in_degree(gr.num_vertices() - 1);

		int32_t lpos = gr.get_vertex_info(0).lpos;
		int32_t rpos = gr.get_vertex_info(gr.num_vertices() - 1).rpos;

		double vv = gr.compute_average_vertex_weight();
		double ee = gr.compute_average_edge_weight();

		printf("subgraph %d, #edges = %lu, #vertices = %lu / %lu, #starting = %d, #ending = %d, range = [%d, %d)\n", 
				i, gr.num_edges(), gr.num_vertices() - 2, root.num_vertices() - 2, d0, dn, lpos, rpos);
	}
	printf("\n");
	return 0;
}
