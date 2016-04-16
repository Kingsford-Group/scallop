#include "splice_graph.h"
#include <sstream>
#include <fstream>
#include <cfloat>
#include <cmath>
#include <algorithm>

using namespace std;

splice_graph::splice_graph()
{}

splice_graph::splice_graph(const splice_graph &gr)
{
	clear();
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		add_vertex();
		set_vertex_string(i, gr.get_vertex_string(i));
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

string splice_graph::get_vertex_string(int v) const
{
	assert(v >= 0 && v < vstr.size());
	return vstr[v];
}

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

int splice_graph::set_vertex_string(int v, string s) 
{
	assert(v >= 0 && v < vv.size());
	if(vstr.size() != vv.size()) vstr.resize(vv.size());
	vstr[v] = s;
	return 0;
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

int splice_graph::clear()
{
	directed_graph::clear();
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

int splice_graph::compute_decomp_paths()
{
	int n = 0;
	for(int i = 0; i < num_vertices(); i++)
	{
		if(degree(i) == 0) continue;
		n++;
	}
	if(n == 0) return 0;
	int m = num_edges();
	return (m - n + 2);
}

int splice_graph::compute_num_paths()
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
			//assert(s < i);
			table[t] += table[s];
		}
	}
	return table[n - 1];
}

bool splice_graph::check_fully_connected()
{
	assert(num_vertices() >= 2);
	if(num_vertices() <= 2) return true;

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

int splice_graph::bfs_w(int s, double w, vector<int> &v, VE &b)
{
	v.assign(num_vertices(), -1);
	b.assign(num_vertices(), null_edge);
	vector<bool> closed;
	closed.resize(num_vertices(), false);
	vector<int> open;
	open.push_back(s);
	v[s] = 0;
	int p = 0;

	while(p < open.size())
	{
		int x = open[p];
		assert(v[x] >= 0);

		p++;
		if(closed[x] == true) continue;
		closed[x] = true;

		edge_iterator it1, it2;
		for(tie(it1, it2) = out_edges(x); it1 != it2; it1++)
		{
			int y = (*it1)->target();
			double ww = get_edge_weight(*it1);
			if(ww < w - SMIN) continue;
			if(v[y] == -1) 
			{
				v[y] = 1 + v[x];
				b[y] = (*it1);
			}
			assert(v[y] <= 1 + v[x]);
			if(closed[y] == true) continue;
			open.push_back(y);
		}
	}
	return 0;
}

int splice_graph::compute_shortest_path_w(int s, int t, double w)
{
	vector<int> v;
	VE b;
	bfs_w(s, w, v, b);
	return v[t];
}

int splice_graph::compute_shortest_path_w(int s, int t, double w, VE &p)
{
	vector<int> v;
	VE b;
	bfs_w(s, w, v, b);
	if(v[t] == -1) return -1;

	p.clear();
	int x = t;
	while(x != s)
	{
		assert(b[x] != null_edge);
		p.push_back(b[x]);
		x = b[x]->source();
	}
	reverse(p.begin(), p.end());
	return v[t];
}

double splice_graph::compute_maximum_path_w(VE &p)
{
	p.clear();
	vector<double> table;		// dynamic programming table
	VE back;					// backtrace edge pointers
	table.resize(num_vertices(), 0);
	back.resize(num_vertices(), null_edge);
	table[0] = DBL_MAX;

	vector<int> tp = topological_sort();
	int n = num_vertices();
	assert(tp.size() == n);
	assert(tp[0] == 0);
	assert(tp[n - 1] == n - 1);

	for(int ii = 1; ii < n; ii++)
	{
		int i = tp[ii];
		if(degree(i) == 0) continue;

		double max_abd = 0;
		edge_descriptor max_edge = null_edge;
		edge_iterator it1, it2;
		for(tie(it1, it2) = in_edges(i); it1 != it2; it1++)
		{
			int s = (*it1)->source();
			int t = (*it1)->target();
			assert(t == i);
			double xw = get_edge_weight(*it1);
			double ww = xw < table[s] ? xw : table[s];
			if(ww >= max_abd)
			{
				max_abd = ww;
				max_edge = *it1;
			}
		}
		assert(max_edge != null_edge);

		back[i] = max_edge;
		table[i] = max_abd;
	}

	int x = n - 1;
	while(true)
	{
		edge_descriptor e = back[x]; 
		if(e == null_edge) break;
		p.push_back(e);
		x = e->source();
	}
	reverse(p.begin(), p.end());

	return table[n - 1];
}

bool splice_graph::compute_optimal_path(VE &p)
{
	p.clear();
	vector<double> table;		// dynamic programming table
	VE back;					// backtrace edge pointers
	table.resize(num_vertices(), 0);
	back.resize(num_vertices(), null_edge);
	table[0] = 0;

	vector<int> tp = topological_sort();
	int n = num_vertices();
	assert(tp.size() == n);
	assert(tp[0] == 0);
	assert(tp[n - 1] == n - 1);

	edge_iterator it1, it2;
	for(int ii = 1; ii < n; ii++)
	{
		int i = tp[ii];
		if(degree(i) == 0) continue;

		double sum = 0;
		for(tie(it1, it2) = in_edges(i); it1 != it2; it1++)
		{
			sum += get_edge_weight(*it1);
		}

		double req = DBL_MAX;
		edge_descriptor ee = null_edge;
		for(tie(it1, it2) = in_edges(i); it1 != it2; it1++)
		{
			int s = (*it1)->source();
			int t = (*it1)->target();
			double xw = get_edge_weight(*it1);
			double ww = sum - xw + table[s];
			if(ww < req)
			{
				req = ww;
				ee = *it1;
			}
		}
		assert(ee != null_edge);

		back[i] = ee;
		table[i] = req;
	}

	edge_descriptor opt = null_edge;
	for(tie(it1, it2) = in_edges(n - 1); it1 != it2; it1++)
	{
		int s = (*it1)->source();
		double w = get_edge_weight(*it1);
		if(w < table[s] + 1) continue;
		opt = *it1;
		break;
	}

	if(opt == null_edge) return false;

	edge_descriptor e = opt;
	while(e != null_edge)
	{
		p.push_back(e);
		int s = e->source();
		e = back[s];
	}
	reverse(p.begin(), p.end());

	return true;
}

double splice_graph::compute_minimum_weight(const VE &p)
{
	double min = DBL_MAX;
	for(int i = 0; i < p.size(); i++)
	{
		double w = get_edge_weight(p[i]);
		if(w < min) min = w;
	}
	return min;
}

int splice_graph::round_weights()
{
	MED m = ewrt;
	for(MED::iterator it = m.begin(); it != m.end(); it++)
	{
		it->second = 0.0;
	}

	while(true)
	{
		VE v;
		double w = compute_maximum_path_w(v);
		if(w <= 0) break;
		double ww = ceil(w);
		
		for(int i = 0; i < v.size(); i++)
		{
			m[v[i]] += ww;
			ewrt[v[i]] -= ww;
			if(ewrt[v[i]] <= 0) ewrt[v[i]] = 0;
		}
	}

	ewrt = m;
	vwrt.assign(num_vertices(), 0);
	edge_iterator it1, it2;
	for(tie(it1, it2) = out_edges(0); it1 != it2; it1++)
	{
		double w = ewrt[*it1];
		vwrt[0] += w;
	}

	for(int i = 1; i < num_vertices(); i++)
	{
		for(tie(it1, it2) = in_edges(i); it1 != it2; it1++)
		{
			double w = ewrt[*it1];
			vwrt[i] += w;
		}
	}

	return 0;
}

int splice_graph::draw(const string &file, const MIS &mis, const MES &mes, double len)
{
	return directed_graph::draw(file, mis, mes, len);
}

int splice_graph::draw(const string &file)
{
	MIS mis;
	char buf[10240];

	for(int i = 0; i < num_vertices(); i++)
	{
		double w = get_vertex_weight(i);
		sprintf(buf, "%.1lf:%s", w, vstr[i].c_str());
		mis.insert(PIS(i, buf));
	}

	MES mes;
	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(); it1 != it2; it1++)
	{
		double w = get_edge_weight(*it1);
		sprintf(buf, "%.2lf", w);
		mes.insert(PES(*it1, buf));
	}
	draw(file, mis, mes, 3.0);
	return 0;
}
