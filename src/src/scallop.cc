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
	round = 0;
	while(iterate());
	return 0;
}

bool scallop::iterate()
{
	while(true)
	{
		print();
		round++;

		int ei;
		vector<int> sub;
		bool flag0 = identify_equation(ei, sub);
		if(flag0 == true) 
		{
			split_edge(ei, sub);
			print();
			round++;
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

				build_adjacent_edges(p);
				print();
				round++;

				connect_adjacent_edges(ex, ey);
				print();
				round++;
			}

			bool b2 = decompose_trivial_vertices();
			if(b2 == true) 
			{
				print();
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

	print();
	round++;

	return true;
}

int scallop::init_super_edges()
{
	mev.clear();
	edge_iterator it1, it2;
	vector<int> v;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
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

			vector<int> v = mev[*it1];
			v.push_back(x);
			v.insert(v.end(), mev[*ot1].begin(), mev[*ot1].end());

			mev.insert(PEV(p, v));
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
		assert(v.size() >= 1);
		xx.push_back(v);
	}
	return xx;
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


int scallop::connect_adjacent_edges(int x, int y)
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
	if(yt == xs) return connect_adjacent_edges(y, x);
	
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
	v.push_back(xt);
	v.insert(v.end(), mev[yy].begin(), mev[yy].end());

	mev.insert(PEV(p, v));

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

int scallop::split_edge(int exi, int eyi)
{
	assert(i2e[exi] != null_edge);
	assert(i2e[eyi] != null_edge);
	edge_descriptor ey = i2e[eyi];
	edge_descriptor ex = i2e[exi];

	double wx = gr.get_edge_weight(ex);
	double wy = gr.get_edge_weight(ey);
	double dx = gr.get_edge_stddev(ex);
	double dy = gr.get_edge_stddev(ey);

	gr.set_edge_weight(ex, wy);
	gr.set_edge_stddev(ex, dx);
	ds.union_set(exi, eyi);

	if(fabs(wx - wy) <= SMIN) return -1;
	assert(wx > wy);

	int s = (ex)->source();
	int t = (ex)->target();

	edge_descriptor p = gr.add_edge(s, t);

	int n = i2e.size();
	i2e.push_back(p);
	assert(e2i.find(p) == e2i.end());
	e2i.insert(PEI(p, n));

	gr.set_edge_weight(p, wx - wy);
	gr.set_edge_stddev(p, dx);
	mev.insert(PEV(p, mev[ex]));

	return n;
}


vector<int> scallop::split_edge(int ei, const vector<int> &sub)
{
	vector<int> v;
	v.push_back(ei);
	int x = ei;
	for(int i = 0; i < sub.size(); i++)
	{
		int y = split_edge(x, sub[i]);
		if(i == sub.size() - 1) assert(y == -1);
		if(y == -1) assert(i == sub.size() - 1);
		if(y == -1) break;
		v.push_back(y);
		x = y;
	}
	return v;
}

bool scallop::identify_equation(int &ei, vector<int> &sub)
{
	vector<int> r = compute_representatives();
	if(r.size() < 2) return false;

	vector<int> x;
	vector<PI> xi;
	for(int i = 0; i < r.size(); i++)
	{
		double w = gr.get_edge_weight(i2e[r[i]]);
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

bool scallop::verify_equation(int ei, const vector<int> &sub)
{
	assert(i2e[ei] != null_edge);
	for(int i = 0; i < sub.size(); i++)
	{
		assert(i2e[sub[i]] != null_edge);
		bool b1 = gr.check_path(i2e[ei], i2e[sub[i]]);
		bool b2 = gr.check_path(i2e[sub[i]], i2e[ei]);
		//printf("check path (%d, %d) = (%c, %c)\n", ei, sub[i], b1 ? 'T' : 'F', b2 ? 'T' : 'F');
		if(b1 == false && b2 == false) return false;
	}
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

int scallop::build_adjacent_edges(const vector<PI> &p)
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

bool scallop::compute_shortest_equal_edges(int &ex, int &ey)
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


int scallop::print()
{
	printf("round %d\n", round);
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
	int n = 0;
	for(int i = 0; i < gr.num_vertices(); i++) 
	{
		if(gr.degree(i) >= 1) n++;
	}
	printf("statistics: %lu edges, %d vertices\n\n", gr.num_edges(), n);

	char buf[1024];
	sprintf(buf, "%s.gr.%d.tex", name.c_str(), round);
	draw_splice_graph(buf);

	sprintf(buf, "%s.nt.%d.tex", name.c_str(), round);
	nt.build(gr);
	nt.draw(buf);

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
