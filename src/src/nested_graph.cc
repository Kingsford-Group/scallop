#include "nested_graph.h"

nested_graph::nested_graph()
{
}

nested_graph::~nested_graph()
{}

int nested_graph::build(directed_graph &gr)
{
	clear();
	mei.clear();
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		add_vertex();
	}

	for(int i = 0; i < gr.num_vertices(); i++)
	{
		if(gr.degree(i) == 0) continue;
		int ip = gr.compute_in_partner(i);
		int op = gr.compute_out_partner(i);
		if(ip == -1 || op == -1) continue;
		edge_descriptor e1 = add_edge(ip, i);
		edge_descriptor e2 = add_edge(i, op);
		edge_descriptor e3 = add_edge(op, ip);
		mei.insert(PEI(e1, i));
		mei.insert(PEI(e2, i));
		mei.insert(PEI(e3, i));
	}
	return 0;
}

int nested_graph::get_in_partner(int x)
{
	edge_iterator it1, it2;
	for(tie(it1, it2) = in_edges(x); it1 != it2; it1++)
	{
		MEI::const_iterator it = mei.find(*it1);
		assert(it != mei.end());
		if(it->second != x) continue;
		return (*it1)->source();
	}
	return -1;
}

int nested_graph::get_out_partner(int x)
{
	edge_iterator it1, it2;
	for(tie(it1, it2) = out_edges(x); it1 != it2; it1++)
	{
		MEI::const_iterator it = mei.find(*it1);
		assert(it != mei.end());
		if(it->second != x) continue;
		return (*it1)->target();
	}
	return -1;
}

vector<int> nested_graph::get_pivots(const vector<int> &p)
{
	vector<int> v;
	if(p.size() <= 1) return v;
	for(int i = 0; i < p.size() - 1; i++)
	{
		vector<edge_descriptor> ve = edges(p[i], p[i + 1]);
		assert(ve.size() >= 1);
		MEI::const_iterator it = mei.find(ve[0]);
		assert(it != mei.end());
		v.push_back(it->second);
	}
	return v;
}

int nested_graph::draw(const string &file) 
{
	MIS mis;
	MES mes;
	for(MEI::iterator it = mei.begin(); it != mei.end(); it++)
	{
		char buf[1024];
		sprintf(buf, "%d", it->second);
		mes.insert(PES(it->first, buf));
	}
	directed_graph::draw(file, mis, mes, 3.0);
	return 0;
}
