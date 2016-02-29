#include "splice_graph.h"
#include "draw.h"

#include <fstream>

int build_splice_graph(const string &file, splice_graph &gr)
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

		add_vertex(gr);
		put(get(vertex_weight, gr), i, weight);
		put(get(vertex_stddev, gr), i, stddev);
	}

	while(fin.getline(line, 10240, '\n'))
	{
		int x, y;
		double weight, stddev;
		stringstream sstr(line);
		sstr>>x>>y>>weight>>stddev;

		assert(x != y);
		assert(x >= 0 && x < num_vertices(gr));
		assert(y >= 0 && y < num_vertices(gr));

		PEB p = add_edge(x, y, gr);
		put(get(edge_weight, gr), p.first, weight);
		put(get(edge_stddev, gr), p.first, stddev);
	}

	fin.close();
	return 0;
}

int draw_splice_graph(const string &file, const splice_graph &gr)
{
	ofstream fout(file.c_str());
	if(fout.fail())
	{
		printf("open file %s error.\n", file.c_str());
		return 0;
	}

	draw_header(fout);

	double len = 1.6;
	fout<<"\\def\\len{"<<len<<"cm}\n";

	// draw vertices
	char sx[1024];
	char sy[1024];
	for(int i = 0; i < num_vertices(gr); i++)
	{
		sprintf(sx, "s%d", i);
		fout.precision(0);
		fout<<fixed;
		fout<<"\\node[mycircle, \\colx, draw, label = below:{";
		//fout<< get(get(vertex_weight, gr), i) << ",";
		fout<< get(get(vertex_weight, gr), i);
		fout<< "}] ("<<sx<<") at ("<<i<<" *\\len, 0.0) {"<<i<<"};\n";
	}

	// draw edges
	adj_iterator ai1, ai2;
	for(int i = 0; i < num_vertices(gr); i++)
	{
		for(tie(ai1, ai2) = adjacent_vertices(i, gr); ai1 != ai2; ai1++)
		{
			int j = *ai1;
			assert(i < j);

			string s;
			char buf[1024];
			out_edge_iterator oi1, oi2;
			for(tie(oi1, oi2) = edge_range(i, j, gr); oi1 != oi2; oi1++)
			{
				double w = get(get(edge_weight, gr), *oi1);
				if(distance(oi1, oi2) == 1) sprintf(buf, "%.0lf", w);
				else sprintf(buf, "%.0lf,", w);

				s.append(buf);
			}

			sprintf(sx, "s%d", i);
			sprintf(sy, "s%d", j);

			double bend = -40;
			if(i + 1 == j) bend = 0;

			fout<<"\\draw[line width = 0.02cm, ->, \\colx, bend right = "<< bend <<"] ("<<sx<<") to node {";
			//fout<< get(get(edge_weight, gr), *it1) <<",";
			fout<< s.c_str() <<"} ("<<sy<<");\n";

		}
	}

	/*
	for(tie(it1, it2) = edges(gr); it1 != it2; it1++)
	{
		int s = source(*it1, gr);
		int t = target(*it1, gr);

		sprintf(sx, "s%d", s);
		sprintf(sy, "s%d", t);
		assert(s < t);
		
		double bend = -40;
		if(s + 1 == t) bend = 0;
		//else if(v[s] % 2 == 0) bend = -30;

		fout<<"\\draw[line width = 0.02cm, ->, \\colx, bend right = "<< bend <<"] ("<<sx<<") to node {";
		//fout<< get(get(edge_weight, gr), *it1) <<",";
		fout<< get(get(edge_weight, gr), *it1) <<"} ("<<sy<<");\n";
	}
	*/

	draw_footer(fout);

	fout.close();
	return 0;
}

int compute_num_paths(const splice_graph &gr)
{
	vector<int> table;
	table.resize(num_vertices(gr), 0);
	table[0] = 1;
	int n = num_vertices(gr);
	for(int i = 1; i < n; i++)
	{
		in_edge_iterator it1, it2;
		for(tie(it1, it2) = in_edges(i, gr); it1 != it2; it1++)
		{
			int s = source(*it1, gr);
			int t = target(*it1, gr);
			assert(t == i);
			assert(s < i);
			table[t] += table[s];
		}
	}
	return table[n - 1];
}

int get_edge_weights(const splice_graph &gr, MED &med)
{
	med.clear();
	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(gr); it1 != it2; it1++)
	{
		double w = get(get(edge_weight, gr), *it1);
		med.insert(PED(*it1, w));
	}
	return 0;
}

int set_edge_weights(splice_graph &gr, const MED &med)
{
	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(gr); it1 != it2; it1++)
	{
		MED::const_iterator it = med.find(*it1);
		put(get(edge_weight, gr), *it1, it->second);
	}
	return 0;
}

int get_vertex_weights(const splice_graph &gr, vector<double> &v)
{
	v.resize(num_vertices(gr), 0);
	for(int i = 0; i < v.size(); i++)
	{
		double w = get(get(vertex_weight, gr), i);
		v[i] = w;
	}
	return 0;
}

int set_vertex_weights(splice_graph &gr, const vector<double> &v)
{
	assert(v.size() == num_vertices(gr));
	for(int i = 0; i < v.size(); i++)
	{
		put(get(vertex_weight, gr), i, v[i]);
	}
	return 0;
}

int get_edge_indices(const splice_graph &gr, VE &i2e, MEI &e2i)
{
	i2e.clear();
	e2i.clear();
	int index = 0;
	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(gr); it1 != it2; it1++)
	{
		e2i.insert(PEI(*it1, index));
		i2e.push_back(*it1);
		index++;
	}
	return 0;
}

bool decide_nested_splice_graph(const splice_graph &gr)
{
	for(int i = 0; i < num_vertices(gr); i++)
	{
		out_edge_iterator it1, it2;
		for(tie(it1, it2) = out_edges(i, gr); it1 != it2; it1++)
		{
			int j = target(*it1, gr);
			assert(j > i);
			for(int k = i + 1; k < j; k++)
			{
				if(exist_direct_path(gr, i, k) == false) continue;
				if(exist_direct_path(gr, k, j) == false) continue;
				out_edge_iterator it3, it4;
				for(tie(it3, it4) = out_edges(k, gr); it3 != it4; it3++)
				{
					int l = target(*it3, gr);
					assert(l > k);
					if(l <= j) continue;
					
					if(exist_direct_path(gr, j, l) == false) continue;
					
					// cross edge found: (i, j) and (k, l)
					bool b1 = unique_path_end(gr, i, k);
					bool b2 = unique_path_start(gr, j, l);
					if(b1 == false && b2 == false) return false;
				}
			}
		}
	}
	return true;
}

bool unique_path_start(const splice_graph &gr, int s, int t)
{
	// assume DAG
	assert(s < t);
	int x = s;
	while(x < t)
	{
		if(out_degree(x, gr) != 1) return false;
		out_edge_iterator it1, it2;
		tie(it1, it2) = out_edges(x, gr);
		int y = target(*it1, gr);
		assert(y > x);
		x = y;
	}
	return true;
}

bool unique_path_end(const splice_graph &gr, int s, int t)
{
	// assume DAG
	assert(s < t);
	int y = t;
	while(y > s)
	{
		if(in_degree(y, gr) != 1) return false;
		in_edge_iterator it1, it2;
		tie(it1, it2) = in_edges(y, gr);
		int x = source(*it1, gr);
		assert(x < y);
		y = x;
	}
	return true;
}

bool exist_direct_path(const splice_graph &gr, int s, int t)
{
	// assume DAG
	assert(s < t);
	vector<bool> closed;
	closed.resize(num_vertices(gr), false);
	vector<int> open;
	open.push_back(s);
	int p = 0;

	while(p < open.size())
	{
		int x = open[p];
		p++;
		if(closed[x] == true) continue;
		closed[x] = true;

		out_edge_iterator it1, it2;
		for(tie(it1, it2) = out_edges(x, gr); it1 != it2; it1++)
		{
			int y = target(*it1, gr);
			if(y == t) return true;
			if(y < t) open.push_back(y);
		}
	}
	return false;
}
