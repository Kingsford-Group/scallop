#include "nested_graph.h"
#include <algorithm>

nested_graph::nested_graph()
{
}

nested_graph::~nested_graph()
{}

int nested_graph::build(directed_graph &gr)
{
	init(gr);
	build_nests(gr);
	build_partners(gr);
	assert(check_nested() == true);
	get_edge_indices(i2e, e2i);
	build_parents();
	build_linkable();
	draw("nested.tex");
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

int nested_graph::build_parents()
{
	parents.assign(num_edges(), -1);
	vector<int> v = topological_sort();
	vector<int> order(num_vertices());
	for(int i = 0; i < v.size(); i++) 
	{
		order[v[i]] = i;
	}
	
	for(int i = 0; i < v.size(); i++)
	{
		int x = v[i];
		if(degree(x) == 0) continue;
		build_parents(x, order);
	}
	return 0;
}

int nested_graph::build_parents(int x, const vector<int> &order)
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

int nested_graph::build_linkable()
{
	linkable.assign(num_edges(), false);
	int n = num_edges();
	for(int i = -1; i < n; i++)
	{
		build_linkable_forward(i);
		build_linkable_backward(i);
	}
	return 0;
}

int nested_graph::build_linkable_forward(int x)
{
	int s = (x == -1) ? 0 : i2e[x]->source();
	set<int> sv;
	vector<int> v;
	edge_iterator it1, it2;
	for(tie(it1, it2) = out_edges(s); it1 != it2; it1++)
	{
		int e = e2i[*it1];
		if(parents[e] != x) continue;
		linkable[e] = true;
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
			if(parents[e] != x) continue;
			if(partners[tt].second != (*it1)->target()) continue;
			assert(linkable[e] == false);
			linkable[e] = true;
			v.push_back(e);
			int t = (*it1)->target();
			assert(sv.find(t) == sv.end());
			sv.insert(s);
		}
	}
	return 0;
}

int nested_graph::build_linkable_backward(int x)
{
	int t = (x == -1) ? num_vertices() - 1 : i2e[x]->target();
	set<int> sv;
	vector<int> v;
	edge_iterator it1, it2;
	for(tie(it1, it2) = in_edges(t); it1 != it2; it1++)
	{
		int e = e2i[*it1];
		if(parents[e] != x) continue;
		if(linkable[e] == true) continue;
		linkable[e] = true;
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
			if(parents[e] != x) continue;
			if(partners[ss].first != (*it1)->source()) continue;
			if(linkable[e] == true) continue;
			linkable[e] = true;
			v.push_back(e);
			int s = (*it1)->source();
			assert(sv.find(s) == sv.end());
			sv.insert(s);
		}
	}
	return 0;
}

bool nested_graph::link(int xs, int xt, int ys, int yt, vector<PI> &p)
{
	if(xt == ys) return true;
	if(yt == xs) return true;
	
	edge_descriptor xe = add_edge(xs, xt);
	edge_descriptor ye = add_edge(ys, yt);

	if(intersect(xe, ye) == true)
	{
		remove_edge(xe);
		remove_edge(ye);
		return false;
	}

	assert(check_nested() == true);
	remove_edge(xe);
	remove_edge(ye);

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
			link(xs, xt, ys, yt, p);
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
		
		if(linkable[e2i[*it1]] == true) ss += "T";
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
		printf("parent of %d = %d\n", i, parents[i]);
	}
	return 0;
}
