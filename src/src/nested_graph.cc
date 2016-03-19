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

int nested_graph::get_in_partner(int x)
{
	return partners[x].first;
}

int nested_graph::get_out_partner(int x)
{
	return partners[x].second;
}

vector<int> nested_graph::get_pivots(const vector<int> &p)
{
	// TODO, might be negative
	vector<int> v;
	return v;
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

		if(partners[s].second == t && partners[s].first >= 0)
		{
			sprintf(buf, "%d", s);
			ss = string(buf);
		}
		if(partners[t].first == s && partners[t].second >= 0)
		{
			sprintf(buf, "%d", t);
			if(ss == "") ss = string(buf);
			else ss += (string(",") + string(buf));
		}
		mes.insert(PES(*it1, ss));
	}
	directed_graph::draw(file, mis, mes, 2.5);
	return 0;
}
