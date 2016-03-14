#include "splice_graph.h"
#include "draw.h"
#include <sstream>
#include <fstream>

using namespace std;

splice_graph::splice_graph()
{}

splice_graph::splice_graph(const splice_graph &gr)
{
	clear();
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		add_vertex();
		set_vertex_weight(i, gr.get_vertex_weight(i));
		set_vertex_stddev(i, gr.get_vertex_stddev(i));
	}

	PEE p = gr.edges();
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		edge_descriptor e = add_edge((*it)->source(), (*it)->target());
		set_edge_weight(e, gr.get_edge_weight(*it));
		set_edge_stddev(e, gr.get_edge_stddev(*it));
	}
}

splice_graph::~splice_graph()
{}

double splice_graph::get_vertex_weight(int v) const
{
	assert(v >= 0 && v < vwrt.size());
	return vwrt[v];
}

double splice_graph::get_vertex_stddev(int v) const
{
	assert(v >= 0 && v < vdev.size());
	return vdev[v];
}

double splice_graph::get_edge_weight(edge_base *e) const
{
	MED::const_iterator it = ewrt.find(e);
	assert(it != ewrt.end());
	return it->second;
}

double splice_graph::get_edge_stddev(edge_base *e) const
{
	MED::const_iterator it = edev.find(e);
	assert(it != edev.end());
	return it->second;
}

int splice_graph::set_vertex_weight(int v, double w) 
{
	assert(v >= 0 && v < vv.size());
	if(vwrt.size() != vv.size()) vwrt.resize(vv.size());
	vwrt[v] = w;
	return 0;
}

int splice_graph::set_vertex_stddev(int v, double w) 
{
	assert(v >= 0 && v < vv.size());
	if(vdev.size() != vv.size()) vdev.resize(vv.size());
	vdev[v] = w;
	return 0;
}

int splice_graph::set_edge_weight(edge_base* e, double w) 
{
	if(ewrt.find(e) != ewrt.end()) ewrt[e] = w;
	else ewrt.insert(PED(e, w));
	return 0;
}

int splice_graph::set_edge_stddev(edge_base* e, double w) 
{
	if(edev.find(e) != edev.end()) edev[e] = w;
	else edev.insert(PED(e, w));
	return 0;
}

MED splice_graph::get_edge_weights() const
{
	return ewrt;
}

vector<double> splice_graph::get_vertex_weights() const
{
	return vwrt;
}

int splice_graph::set_edge_weights(const MED &med)
{
	ewrt = med;
	return 0;
}

int splice_graph::set_vertex_weights(const vector<double> &v)
{
	vwrt = v;
	return 0;
}

int splice_graph::get_edge_indices(VE &i2e, MEI &e2i) const
{
	i2e.clear();
	e2i.clear();
	int index = 0;
	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(); it1 != it2; it1++)
	{
		e2i.insert(PEI(*it1, index));
		i2e.push_back(*it1);
		index++;
	}
	return 0;
}

int splice_graph::clear()
{
	graph_base::clear();
	vwrt.clear();
	vdev.clear();
	ewrt.clear();
	edev.clear();
	return 0;
}

int splice_graph::build(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail()) 
	{
		printf("open file %s error\n", file.c_str());
		return 0;
	}

	char line[10240];
	// get the number of vertices
	fin.getline(line, 10240, '\n');	
	int n = atoi(line);

	for(int i = 0; i < n; i++)
	{
		double weight, stddev;
		fin.getline(line, 10240, '\n');	
		stringstream sstr(line);
		sstr>>weight>>stddev;

		add_vertex();
		set_vertex_weight(i, weight);
		set_vertex_stddev(i, stddev);
	}

	while(fin.getline(line, 10240, '\n'))
	{
		int x, y;
		double weight, stddev;
		stringstream sstr(line);
		sstr>>x>>y>>weight>>stddev;

		assert(x != y);
		assert(x >= 0 && x < num_vertices());
		assert(y >= 0 && y < num_vertices());

		edge_descriptor p = add_edge(x, y);
		set_edge_weight(p, weight);
		set_edge_stddev(p, stddev);
	}

	fin.close();
	return 0;
}

int splice_graph::write(const string &file) const
{
	ofstream fin(file.c_str());
	if(fin.fail()) 
	{
		printf("open file %s error\n", file.c_str());
		return 0;
	}
	
	fin<<fixed;
	fin.precision(2);
	int n = num_vertices();
	
	fin<<n<<endl;
	for(int i = 0; i < n; i++)
	{
		double weight = get_vertex_weight(i);
		double stddev = get_vertex_stddev(i);
		fin<<weight<<" "<<stddev<<endl;
	}

	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(); it1 != it2; it1++)
	{
		int s = (*it1)->source(); 
		int t = (*it1)->target();
		double weight = get_edge_weight(*it1);
		double stddev = get_edge_stddev(*it1);
		fin<<s<<" "<<t<<" "<<weight<<" "<<stddev<<endl;
	}
	fin.close();
	return 0;
}

int splice_graph::simulate(int n, int m)
{
	clear();
	for(int i = 0; i < n; i++)
	{
		add_vertex();
		set_vertex_weight(i, 1);
		set_vertex_stddev(i, 1);
	}

	for(int i = 0; i < m; i++)
	{
		int s = rand() % (n - 1);
		int t = s + 1 + rand() % (n - s - 1);
		edge_descriptor p = add_edge(s, t);
		set_edge_weight(p, 1);
		set_edge_stddev(p, 1);
	}
	return 0;
}

int splice_graph::compute_num_paths() const
{
	vector<int> table;
	int n = num_vertices();
	table.resize(n, 0);
	table[0] = 1;
	for(int i = 1; i < n; i++)
	{
		edge_iterator it1, it2;
		for(tie(it1, it2) = in_edges(i); it1 != it2; it1++)
		{
			int s = (*it1)->source();
			int t = (*it1)->target();
			assert(t == i);
			assert(s < i);
			table[t] += table[s];
		}
	}
	return table[n - 1];
}

bool splice_graph::check_nested() const
{
	vector< vector<int> > vv;
	vv.resize(num_vertices());
	for(int i = 0; i < vv.size(); i++) bfs(i, vv[i]);

	for(int i = 0; i < num_vertices(); i++)
	{
		edge_iterator it1, it2;
		for(tie(it1, it2) = out_edges(i); it1 != it2; it1++)
		{
			int j = (*it1)->target();
			assert(j > i);
			for(int k = i + 1; k < j; k++)
			{
				if(vv[i][k] == -1) continue;
				if(vv[k][j] == -1) continue;
				edge_iterator it3, it4;
				for(tie(it3, it4) = out_edges(k); it3 != it4; it3++)
				{
					int l = (*it3)->target();
					assert(l > k);
					if(l <= j) continue;
					
					if(vv[j][l] == -1) continue;
					return false;
				}
			}
		}
	}
	return true;
}

bool splice_graph::check_fully_connected() const
{
	if(num_vertices() <= 1) return true;

	vector<int> s;
	vector<int> t;
	bfs(0, s);
	bfs_reverse(num_vertices() - 1, t);

	for(int i = 0; i < s.size(); i++)
	{
		if(s[i] == -1) return false;
	}

	for(int i = 0; i < t.size(); i++)
	{
		if(t[i] == -1) return false;
	}

	return true;
}
