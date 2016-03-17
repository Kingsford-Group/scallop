#include "scallop3.h"
#include "subsetsum.h"
#include <cstdio>
#include <iostream>
#include <cfloat>

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
	round = 0;
	while(iterate());
	return 0;
}

bool scallop3::iterate()
{
	char buf[1024];
	while(true)
	{
		printf("round %d, start\n", round);
		print();

		sprintf(buf, "%s.gr.%d.tex", name.c_str(), round);
		this->draw_splice_graph(buf);
		round++;

		int ei;
		vector<int> sub;
		bool flag0 = identify_equation(ei, sub);
		if(flag0 == true) 
		{
			split_edge(ei, sub);

			printf("round %d, split edge %d\n", round, ei);
			print();

			sprintf(buf, "%s.gr.%d.tex", name.c_str(), round);
			this->draw_splice_graph(buf);
			round++;
		}

		bool flag1 = false;
		while(true)
		{
			compute_intersecting_edges();

			int ex, ey;
			vector<int> p;
			bool b1 = identify_linkable_edges(ex, ey, p);
			if(b1 == true)
			{
				assert(p.size() >= 1);
				printf("linkable edges = (%d, %d), path = (%d", ex, ey, p[0]);
				for(int i = 1; i < p.size(); i++) printf(", %d", p[i]);
				printf(")\n");
				assert(ex >= 0 && ey >= 0);

				build_adjacent_edges(ex, ey, p);
				connect_adjacent_edges(ex, ey);

				printf("round %d, connect edge %d and %d\n", round, ex, ey);
				print();

				sprintf(buf, "%s.gr.%d.tex", name.c_str(), round);
				this->draw_splice_graph(buf);
				round++;
			}

			bool b2 = decompose_trivial_vertices();
			if(b2 == true) 
			{
				printf("round %d, decompose trivial vertex\n", round);
				print();

				sprintf(buf, "%s.gr.%d.tex", name.c_str(), round);
				this->draw_splice_graph(buf);
				round++;
			}

			if(b1 == true || b2 == true) flag1 = true;
			if(b1 == false && b2 == false) break;
		}

		if(flag0 == false && flag1 == false) break;
	}

	int ex, ey;
	bool b = compute_shortest_equal_edges(ex, ey);

	if(b == false) return false;

	printf("shortest equal path (%d, %d)\n", ex, ey);

	connect_equal_edges(ex, ey);

	printf("round %d, connect %d and %d\n", round, ex, ey);
	print();

	sprintf(buf, "%s.gr.%d.tex", name.c_str(), round);
	this->draw_splice_graph(buf);
	round++;

	return true;
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
			bool b = init_trivial_vertex(i);
			if(b == true) flag = true;
		}
		if(flag == false) break;
	}
	return 0;
}

bool scallop3::init_trivial_vertex(int x)
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

int scallop3::connect_equal_edges(int x, int y)
{
	assert(i2e[x] != null_edge);
	assert(i2e[y] != null_edge);

	edge_descriptor xx = i2e[x];
	edge_descriptor yy = i2e[y];

	int xs = source(xx, gr);
	int xt = target(xx, gr);
	int ys = source(yy, gr);
	int yt = target(yy, gr);

	vector<int> p;
	bool b = gr.compute_shortest_path(xt, ys, p);
	assert(b == true);
	assert(p.size() >= 1);

	double wx = gr.get_edge_weight(xx);
	double wy = gr.get_edge_weight(yy);
	assert(fabs(wx - wy) <= SMIN);

	double wp = gr.compute_bottleneck_weight(p);
	assert(wp >= wx - SMIN);

	int pre = x;
	for(int i = 0; i < p.size() - 1; i++)
	{
		edge_iterator it1, it2;
		int ee = -1;
		double dist = DBL_MAX;
		for(tie(it1, it2) = gr.out_edges(p[i]); it1 != it2; it1++)
		{
			if((*it1)->target() != p[i + 1]) continue;
			double w = gr.get_edge_weight(*it1);
			if(w <= wx - SMIN) continue;
			double d = fabs(w - wx);
			if(d < dist)
			{
				dist = d;
				ee = e2i[*it1];
			}
		}
		assert(ee != -1);

		split_edge(ee, pre);
		int e1 = connect_adjacent_edges(pre, ee);
		assert(e1 != -1);
		pre = e1;
	}

	split_edge(y, pre);
	connect_adjacent_edges(pre, y);

	return 0;
}


int scallop3::connect_adjacent_edges(int x, int y)
{
	if(i2e[x] == null_edge) return -1;
	if(i2e[y] == null_edge) return -1;

	edge_descriptor xx = i2e[x];
	edge_descriptor yy = i2e[y];

	int xs = source(xx, gr);
	int xt = target(xx, gr);
	int ys = source(yy, gr);
	int yt = target(yy, gr);

	if(xt != ys && yt != xs) return -1;
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

	return n;
}

int scallop3::split_edge(int exi, int eyi)
{
	assert(i2e[exi] != null_edge);
	assert(i2e[eyi] != null_edge);
	edge_descriptor ey = i2e[eyi];
	edge_descriptor ex = i2e[exi];

	double wx = get_edge_weight(ex, gr);
	double wy = get_edge_weight(ey, gr);
	double dx = get_edge_stddev(ex, gr);
	double dy = get_edge_stddev(ey, gr);

	set_edge_weight(ex, wy, gr);
	set_edge_stddev(ex, dx, gr);
	ds.union_set(exi, eyi);

	if(fabs(wx - wy) <= SMIN) return -1;
	assert(wx > wy);

	int s = source(ex, gr);
	int t = target(ex, gr);

	PEB p = add_edge(s, t, gr);
	assert(p.second == true);

	int n = i2e.size();
	i2e.push_back(p.first);
	assert(e2i.find(p.first) == e2i.end());
	e2i.insert(PEI(p.first, n));

	set_edge_weight(p.first, wx - wy, gr);
	set_edge_stddev(p.first, dx, gr);
	mev.insert(PEV(p.first, mev[ex]));

	return n;
}


vector<int> scallop3::split_edge(int ei, const vector<int> &sub)
{
	vector<int> v;
	v.push_back(ei);
	int x = ei;
	for(int i = 0; i < sub.size(); i++)
	{
		printf("split edge %d w.r.t. %d\n", x, sub[i]);
		int y = split_edge(x, sub[i]);
		if(i == sub.size() - 1) assert(y == -1);
		if(y == -1) assert(i == sub.size() - 1);
		if(y == -1) break;
		v.push_back(y);
		x = y;
	}
	return v;
}

bool scallop3::identify_equation(int &ei, vector<int> &sub)
{
	vector<int> r = compute_representatives();
	if(r.size() < 2) return false;

	vector<int> x;
	vector<PI> xi;
	for(int i = 0; i < r.size(); i++)
	{
		double w = get_edge_weight(i2e[r[i]], gr);
		x.push_back((int)(w));
		xi.push_back(PI((int)(w), i));
	}

	sort(xi.begin(), xi.end());

	int min = -1;
	for(int i = 1; i < xi.size(); i++)
	{
		vector<int> v;
		for(int j = 0; j <= i; j++) 
		{
			v.push_back(xi[j].first);
		}

		subsetsum sss(v);
		sss.solve();

		int eix = r[xi[i].second];
		vector<int> ssx;
		for(int j = 0; j < sss.subset.size(); j++)
		{
			int k = sss.subset[j];
			int kk = xi[k].second;
			ssx.push_back(r[kk]);
		}

		bool b = verify_equation(eix, ssx);
		if(b == false) continue;

		int err = (int)fabs(sss.opt - xi[i].first);
		if(min == -1 || err < min)
		{
			min = err;
			ei = eix;
			sub = ssx;
		}
		else if(err == min && sub.size() > ssx.size())
		{
			min = err;
			ei = eix;
			sub = ssx;
		}
	}

	if(min == -1) return false;

	printf("%s closest subset for edge %d:%.0lf has %lu edges, error = %d, subset = (", name.c_str(), ei, gr.get_edge_weight(i2e[ei]), sub.size(), min);
	for(int i = 0; i < sub.size() - 1; i++) printf("%d:%.0lf, ", sub[i], gr.get_edge_weight(i2e[sub[i]]));
	printf("%d:%.0lf)\n", sub[sub.size() - 1], gr.get_edge_weight(i2e[sub[sub.size() - 1]]));
	
	if(min >= 1) return false;
	else return true;
}

bool scallop3::verify_equation(int ei, const vector<int> &sub)
{
	assert(i2e[ei] != null_edge);
	for(int i = 0; i < sub.size(); i++)
	{
		assert(i2e[sub[i]] != null_edge);
		bool b1 = gr.check_directed_path(i2e[ei], i2e[sub[i]]);
		bool b2 = gr.check_directed_path(i2e[sub[i]], i2e[ei]);
		//printf("check path (%d, %d) = (%c, %c)\n", ei, sub[i], b1 ? 'T' : 'F', b2 ? 'T' : 'F');
		if(b1 == false && b2 == false) return false;
	}
	return true;
}

int scallop3::compute_intersecting_edges()
{
	sis.clear();
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		for(int j = i + 1; j < i2e.size(); j++)
		{
			if(i2e[j] == null_edge) continue;
			bool b = gr.intersect(i2e[i], i2e[j]);
			if(b == false) continue;
			int s1 = i2e[i]->source();
			int t1 = i2e[i]->target();
			int s2 = i2e[j]->source();
			int t2 = i2e[j]->target();
			sis.insert(PI(s1, t1));
			sis.insert(PI(s2, t2));
		}
	}
	return 0;
}

bool scallop3::check_linkable(int ex, int ey, vector<int> &p)
{
	assert(i2e[ex] != null_edge);
	assert(i2e[ey] != null_edge);
	bool b1 = gr.check_directed_path(i2e[ex], i2e[ey]);
	bool b2 = gr.check_directed_path(i2e[ey], i2e[ex]);
	assert(b1 == false || b2 == false);
	if(b1 == false && b2 == false) return false;
	if(b2 == true) return check_linkable(ey, ex, p);

	bool b = gr.compute_shortest_path(i2e[ex], i2e[ey], p);
	if(b == false) return false;
	assert(p.size() >= 1);
	if(p.size() == 1) return true;
	
	for(int i = 0; i < p.size() - 1; i++)
	{
		PI pi(p[i], p[i + 1]);
		if(sis.find(pi) != sis.end()) return false;
	}

	int li = 0;
	int ri = p.size() - 1;
	while(li < ri)
	{
		int l1 = p[li];
		int r1 = p[ri];
		int l2 = p[li + 1];
		int r2 = p[ri - 1];

		int lr = gr.compute_out_ancestor(l1);
		int rr = gr.compute_out_ancestor(r1);
		int ll = gr.compute_in_ancestor(l1);
		int rl = gr.compute_in_ancestor(r1);

		if(lr == l2 && sis.find(PI(ll, l1)) == sis.end())
		{
			p[li] = 0 - p[li];
			li++;
		}
		else if(rl == r2 && sis.find(PI(r1, rr)) == sis.end())
		{
			p[ri] = 0 - p[ri];
			ri--;
		}
		else return false;
	}
	return true;
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
				bool b = check_linkable(v[j], v[k], p);
				if(b == true)
				{
					bool b1 = gr.check_directed_path(i2e[v[j]], i2e[v[k]]);
					bool b2 = gr.check_directed_path(i2e[v[k]], i2e[v[j]]);
					assert(b1 != b2);
					if(b1 == true)
					{
						ex = v[j];
						ey = v[k];
					}
					else
					{
						ex = v[k];
						ey = v[j];
					}
					flag = true;
					break;
				}
			}
			if(flag == true) break;
		}
		if(flag == true) break;
	}
	if(flag == false) return false;
	return true;
}

int scallop3::build_adjacent_edges(int ex, int ey, const vector<int> &p)
{
	int li = 0;
	int ri = p.size() - 1;
	while(li < ri)
	{
		int l1 = (int)fabs(p[li]);
		int r1 = (int)fabs(p[ri]);
		int l2 = (int)fabs(p[li + 1]);
		int r2 = (int)fabs(p[ri - 1]);

		if(p[li] < 0)
		{
			int lr = gr.compute_out_ancestor(l1);
			int ll = gr.compute_in_ancestor(l1);
			assert(lr == l2);

			gr.exchange(ll, l1, l2);

			printf("exchange %d %d %d\n", ll, l1, l2);
			char buf[1024];
			sprintf(buf, "%d.%d.%d.tex", ll, l1, l2);
			draw_splice_graph(buf);

			li++;
		}
		else if(p[ri] < 0)
		{
			int rl = gr.compute_in_ancestor(r1);
			int rr = gr.compute_out_ancestor(r1);
			assert(rl == r2);

			gr.exchange(r2, r1, rr);

			printf("exchange %d %d %d\n", r2, r1, rr);
			char buf[1024];
			sprintf(buf, "%d.%d.%d.tex", r2, r1, rr);
			draw_splice_graph(buf);

			ri--;
		}
		else assert(false);
	}
	return 0;
}

bool scallop3::decompose_trivial_vertices()
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
				connect_adjacent_edges(v[k], sub[k]);
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
				connect_adjacent_edges(sub[k], v[k]);
			}
			flag = true;
		}
	}
	return flag;
}

bool scallop3::compute_shortest_equal_edges(int &ex, int &ey)
{
	ex = ey = -1;
	vector< vector<int> > vv = compute_disjoint_sets();
	bool flag = false;
	int min = -1;

	for(int i = 0; i < vv.size(); i++)
	{
		vector<int> &v = vv[i];
		if(v.size() == 1) continue;
		for(int j = 0; j < v.size(); j++)
		{
			double wx = gr.get_edge_weight(i2e[v[j]]);
			for(int k = j + 1; k < v.size(); k++)
			{
				vector<int> p1;
				vector<int> p2;
				bool b1 = gr.compute_shortest_path(i2e[v[j]], i2e[v[k]], p1);
				bool b2 = gr.compute_shortest_path(i2e[v[k]], i2e[v[j]], p2);

				assert(b1 == false || b2 == false);
				if(b1 == false && b2 == false) continue;
				assert(b1 != b2);

				if(b1 == true)
				{
					double w = gr.compute_bottleneck_weight(p1);
					if(wx - SMIN >= w) continue;
					if(min == -1 || min > p1.size())
					{
						min = p1.size();
						ex = v[j];
						ey = v[k];
					}
				}
				else
				{
					double w = gr.compute_bottleneck_weight(p2);
					if(wx - SMIN >= w) continue;
					if(min == -1 || min > p2.size())
					{
						min = p2.size();
						ex = v[k];
						ey = v[j];
					}
				}
			}
		}
	}
	if(min == -1) return false;
	else return true;
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

		if(v.size() == 1) continue;
		int w = (int)(get_edge_weight(i2e[v[0]], gr));

		printf("edge set %d, weight = %d, #edges = %lu, set = (%d", i, w, v.size(), v[0]);
		for(int j = 1; j < v.size(); j++) printf(", %d", v[j]);
		printf(")\n");
	}
	int n = 0;
	for(int i = 0; i < gr.num_vertices(); i++) 
	{
		if(gr.degree(i) >= 1) n++;
	}
	printf("statistics: %lu edges, %d vertices\n", gr.num_edges(), n);

	// print in and out partners
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		if(gr.degree(i) == 0) continue;
		int ip = gr.compute_in_partner(i);
		int op = gr.compute_out_partner(i);
		printf("partner of %d: in = %d, out = %d\n", i, ip, op);
	}
	printf("\n");
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

	gr.draw(file, mis, mes, 4.0);
	return 0;
}

