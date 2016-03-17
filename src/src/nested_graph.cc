#include "nested_graph.h"

nested_graph::nested_graph(const graph_base &gr)
{
	build(gr);
}

nested_graph::~nested_graph()
{}

int nested_graph::build(const graph_base &gr)
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
	graph_base::draw(file, mis, mes, 3.0);
	return 0;
}
