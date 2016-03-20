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
	//test_linking();
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

int nested_graph::bfs_search(int s, vector<int> &table, vector<int> &open)
{
	table.assign(num_vertices() * 2, -1);
	
	open.clear();
	open.push_back(s);
	table[s] = 0;
	int p = 0;

	while(p < open.size())
	{
		int x = open[p];
		assert(table[x] >= 0);
		p++;

		edge_iterator it1, it2;
		if(x < num_vertices())
		{
			for(tie(it1, it2) = in_edges(x); it1 != it2; it1++)
			{
				int ss = (*it1)->source();
				int tt = (*it1)->target();
				assert(tt == x);

				int y = ss + num_vertices();
				if(table[y] >= 0) continue;

				assert(x != y);
				table[y] = x;
				open.push_back(y);
			}
		}
		else if(x < 2 * num_vertices())
		{
			for(tie(it1, it2) = out_edges(x - num_vertices()); it1 != it2; it1++)
			{
				int ss = (*it1)->source();
				int tt = (*it1)->target();
				assert(ss == x - num_vertices());

				int y = tt;
				if(table[y] >= 0) continue;

				assert(x != y);
				table[y] = x;
				open.push_back(y);
			}
		}
		else assert(false);

		if(partners[x].first == -1 || partners[x].second == -1) continue;

		if(x < num_vertices())
		{
			for(tie(it1, it2) = out_edges(x); it1 != it2; it1++)
			{
				int ss = (*it1)->source();
				int tt = (*it1)->target();
				assert(ss == x);

				if(partners[x].second != tt) continue;

				int y = tt;
				if(table[y] >= 0) continue;

				assert(x != y);
				table[y] = x;
				open.push_back(y);
			}
		}
		else if(x < 2 * num_vertices())
		{
			for(tie(it1, it2) = in_edges(x - num_vertices()); it1 != it2; it1++)
			{
				int ss = (*it1)->source();
				int tt = (*it1)->target();
				assert(tt == x - num_vertices());

				if(partners[x].first != ss) continue;

				int y = ss + num_vertices();
				if(table[y] >= 0) continue;

				assert(x != y);
				table[y] = x;
				open.push_back(y);
			}
		}
		else assert(false);
	}
	table[s] = -1;
	return 0;
}

bool nested_graph::link(int s, int t, vector<PI> &p)
{
	int n = num_vertices();

	vector<int> s1, s2;
	vector<int> t1, t2;
	bfs_search(s, s1, s2);
	bfs_search(t, t1, t2);

	int ks = -1, kt = -1;
	for(int i = 0; i < t2.size(); i++)
	{
		int x = t2[i];
		int y = (x >= n ? x - n : x + n);
		if(s1[y] == -1) continue;
		ks = y;
		kt = x;
		break;
	}

	if(ks == -1 || kt == -1) return false;

	p.clear();
	
	int x = ks;
	while(true)
	{
		int y = s1[x];
		if(y == -1) break;
		if(x < n && y < n) p.push_back(PI(y, -1));
		else if(x >= n && y >= n) p.push_back(PI(y - n, -1));
		else if(x >= n && y < n) p.push_back(PI(x - n, y));
		else p.push_back(PI(y - n, x));
		x = y;
	}

	x = kt;
	while(true)
	{
		int y = t1[x];
		if(y == -1) break;
		if(x < n && y < n) p.push_back(PI(y, -1));
		else if(x >= n && y >= n) p.push_back(PI(y - n, -1));
		else if(x >= n && y < n) p.push_back(PI(x - n, y));
		else p.push_back(PI(y - n, x));
		x = y;
	}

	reverse(p.begin(), p.end());
	return true;
}

int nested_graph::test_linking()
{
	for(int i = 0; i < num_vertices() * 2; i++)
	{
		for(int j = 0; j < num_vertices() * 2; j++)
		{
			vector<PI> p;
			bool b = link(i, j, p);
			if(b == false) continue;

			printf(" connect (%d, %d): ", i, j);
			for(int k = 0; k < p.size(); k++) printf("(%d, %d) ", p[k].first, p[k].second);
			printf("\n");
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
