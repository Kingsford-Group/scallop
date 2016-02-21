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
	edge_iterator it1, it2;
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

