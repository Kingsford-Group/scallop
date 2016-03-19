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
	build_minimal_intervals(gr);
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

int nested_graph::build_minimal_intervals(directed_graph &gr)
{
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		if(gr.degree(i) == 0) continue;
		for(int j = 0; j < gr.num_vertices(); j++)
		{
			if(gr.degree(j) == 0) continue;
			int b = gr.check_partner(i, j);
			if(b == -1) continue;
			add_edge(i, j);
		}
	}
	return 0;
}

int nested_graph::build_all_intervals(directed_graph &gr)
{
	vector< set<int> > vvi(gr.num_vertices());
	vector< set<int> > vvo(gr.num_vertices());
	vector< set<edge_descriptor> > vei(gr.num_vertices());
	vector< set<edge_descriptor> > veo(gr.num_vertices());
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		gr.compute_in_content(i, vvi[i], vei[i]);
		gr.compute_out_content(i, vvo[i], veo[i]);
	}

	vector<int> vv(gr.num_vertices());
	vector<edge_descriptor> ve(gr.num_edges());
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		for(int j = 0; j < gr.num_vertices(); j++)
		{
			if(i == j) continue;
			vector<int>::iterator iv = set_intersection(vvo[i].begin(), vvo[i].end(), vvi[j].begin(), vvi[j].end(), vv.begin());
			vector<edge_descriptor>::iterator ie = set_intersection(veo[i].begin(), veo[i].end(), vei[j].begin(), vei[j].end(), ve.begin());
			vector<int> v(vv.begin(), iv);
			set<edge_descriptor> se(ve.begin(), ie);

			printf("checking %d and %d: %lu internal vertices, %lu internal edges\n", i, j, v.size(), se.size());
			if(se.size() == 0) continue;

			bool b = true;
			for(int k = 0; k < v.size(); k++)
			{
				printf(" internal vertex %d\n", v[k]);
				edge_iterator it1, it2;
				for(tie(it1, it2) = gr.in_edges(v[k]); it1 != it2; it1++)
				{
					if(se.find(*it1) == se.end()) b = false;
					if(b == false) break;
				}
				if(b == false) break;
				for(tie(it1, it2) = gr.out_edges(v[k]); it1 != it2; it1++)
				{
					if(se.find(*it1) == se.end()) b = false;
					if(b == false) break;
				}
				if(b == false) break;
			}
			if(b == false) continue;

			edge_descriptor e = add_edge(i, j);
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
