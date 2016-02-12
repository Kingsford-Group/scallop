#include "sgraph.h"
#include "draw.h"

#include <iomanip>
#include <cfloat>

sgraph::sgraph(const bbase &b)
	:bundle(b)
{}

sgraph::~sgraph()
{}

int sgraph::solve()
{
	build_graph();
	check();
	build_paths();
	return 0;
}

int sgraph::build_graph()
{
	// vertices: start, each region, end
	add_vertex(gr);
	put(get(vertex_weight, gr), 0, 0);
	for(int i = 0; i < regions.size(); i++)
	{
		add_vertex(gr);
		put(get(vertex_weight, gr), i + 1, regions[i].ave_abd);
	}
	add_vertex(gr);
	put(get(vertex_weight, gr), regions.size() + 1, 0);

	// edges: connecting adjacent regions => e2w
	for(int i = 0; i < regions.size() - 1; i++)
	{
		region &x = regions[i];
		region &y = regions[i + 1];

		if(x.empty || y.empty) continue;

		if(x.right_break()) continue;
		if(y.left_break()) continue;
		if(x.rtype == RIGHT_BOUNDARY) continue;
		if(y.ltype == LEFT_BOUNDARY) continue;

		assert(x.rpos == y.lpos);
		int32_t xr = compute_overlap(imap, x.rpos - 1);
		int32_t yl = compute_overlap(imap, y.lpos);

		double w = xr < yl ? xr : yl;
		//double w = x.ave_abd < y.ave_abd ? x.ave_abd : y.ave_abd;

		PEB p = add_edge(i + 1, i + 2, gr);
		assert(p.second == true);
		put(get(edge_weight, gr), p.first, w);
	}

	// edges: each bridge => and e2w
	for(int i = 0; i < bridges.size(); i++)
	{
		bridge &b = bridges[i];
		PEB p = add_edge(b.lrgn + 1, b.rrgn + 1, gr);
		assert(p.second == true);
		put(get(edge_weight, gr), p.first, b.count);

		//assert(e2b.find(p.first) == e2b.end());
		//e2b.insert(PEI(p.first, i));
	}

	// edges: connecting start/end and regions
	int ss = 0;
	int tt = regions.size() + 1;
	for(int i = 0; i < regions.size(); i++)
	{
		region &r = regions[i];
		if(r.empty == true) continue;

		// TODO
		//if(r.ltype == LEFT_BOUNDARY || r.ltype == START_BOUNDARY)
		if(r.left_break() || r.ltype == LEFT_BOUNDARY || r.ltype == START_BOUNDARY)
		{
			PEB p = add_edge(ss, i + 1, gr);
			assert(p.second == true);
			put(get(edge_weight, gr), p.first, r.ave_abd);
		}

		// TODO
		//if(r.rtype == RIGHT_BOUNDARY || r.rtype == END_BOUNDARY) 
		if(r.right_break() || r.rtype == RIGHT_BOUNDARY || r.rtype == END_BOUNDARY) 
		{
			PEB p = add_edge(i + 1, tt, gr);
			assert(p.second == true);
			put(get(edge_weight, gr), p.first, r.ave_abd);
		}
	}

	return 0;
}

int sgraph::build_paths()
{
	path p;
	compute_maximum_path(p);
	paths.push_back(p);
	return 0;
}

int sgraph::compute_maximum_path(path &p)
{
	// TODO: use the weight on vertices right now
	vector<double> table;		// dynamic programming table
	vector<int> back;			// back pointers
	table.resize(num_vertices(gr), 0);
	back.resize(num_vertices(gr), -1);
	table[0] = DBL_MAX;
	int n = num_vertices(gr);
	for(int i = 1; i < n; i++)
	{
		double abd = get(get(vertex_weight, gr), i);
		if(i == n - 1) abd = DBL_MAX;

		double max_abd = 0;
		int max_idx = -1;
		in_edge_iterator it1, it2;
		for(tie(it1, it2) = in_edges(i, gr); it1 != it2; it1++)
		{
			int s = source(*it1, gr);
			int t = target(*it1, gr);
			assert(t == i);
			assert(s < i);
			if(table[s] >= max_abd)
			{
				max_abd = table[s];
				max_idx = s;
			}
		}

		back[i] = max_idx;
		table[i] = max_abd < abd ? max_abd : abd;
	}

	if(table[n - 1] <= 0) return 0;

	p.abd = table[n - 1];
	p.v.clear();

	int b = n - 1;
	while(true)
	{
		p.v.push_back(b);
		if(b == 0) break;
		b = back[b];
		assert(b != -1);
	}
	reverse(p.v.begin(), p.v.end());

	return 0;
}

PEB sgraph::get_max_in_edge(int x)
{
	in_edge_iterator it1, it2;
	double max = -1;
	PEB p;
	p.second = false;
	for(tie(it1, it2) = in_edges(x, gr); it1 != it2; it1++)
	{
		double w =get(get(edge_weight, gr), *it1);
		if(w > max)
		{
			p.first = *it1;
			p.second = true;
			max = w;
		}
	}
	return p;
}

PEB sgraph::get_max_out_edge(int x)
{
	out_edge_iterator it1, it2;
	double max = -1;
	PEB p;
	p.second = false;
	for(tie(it1, it2) = out_edges(x, gr); it1 != it2; it1++)
	{
		double w = get(get(edge_weight, gr), *it1);
		if(w > max)
		{
			p.first = *it1;
			p.second = true;
			max = w;
		}
	}
	return p;
}

int sgraph::check()
{
	// check
	if(regions.size() == 1) return 0;
	for(int i = 0; i < regions.size(); i++)
	{
		region &r = regions[i];
		if(r.empty == false) continue;
		assert(r.empty == true);

		if(in_degree(i + 1, gr) != 0) printf("problem at region %d\n", i);
		assert(in_degree(i + 1, gr) == 0);

		if(out_degree(i + 1, gr) != 0) printf("problem at region %d\n", i);
		assert(out_degree(i + 1, gr) == 0);
	}

	return 0;
}

int sgraph::print(int index)
{
	printf("Bundle %d: ", index);
	printf("tid = %d, #hits = %lu, range = %s:%d-%d\n", tid, hits.size(), chrm.c_str(), lpos, rpos);
	// print hits
	/*
	for(int i = 0; i < hits.size(); i++)
	{
		hits[i].print();
	}
	*/

	// print bridges 
	for(int i = 0; i < bridges.size(); i++)
	{
		bridges[i].print(i);
	}

	// print boundaries
	for(int i = 0; i < boundaries.size(); i++)
	{
		boundaries[i].print(i);
	}

	// print regions
	for(int i = 0; i < regions.size(); i++)
	{
		regions[i].print(i);
	}

	// print paths
	for(int i = 0; i < paths.size(); i++)
	{
		paths[i].print(i);
	}

	printf("\n");
	return 0;
}

int sgraph::draw(const string &file)
{
	ofstream fout(file.c_str());
	if(fout.fail())
	{
		printf("open file %s error.\n", file.c_str());
		return 0;
	}

	draw_header(fout);

	double len = 1.5;
	fout<<"\\def\\len{"<<len<<"cm}\n";

	vector<int> v;
	int vi = 1;
	v.push_back(0);
	for(int i = 0; i < regions.size(); i++)
	{
		v.push_back(vi++);
		//if(regions[i].empty == true) v.push_back(-1);
		//else v.push_back(vi++);
	}
	v.push_back(vi);

	assert(v.size() == num_vertices(gr));

	// draw vertices
	char sx[1024];
	char sy[1024];
	for(int i = 0; i < num_vertices(gr); i++)
	{
		if(v[i] == -1) continue;
		sprintf(sx, "s%d", i);
		double px = v[i] * len;
		double py = 0.0;
		fout.precision(1);
		fout<<fixed;
		fout<<"\\node[mycircle, \\colx, draw, label = below:{"<< get(get(vertex_weight, gr), i) << "}] ("<<sx<<") at ("<<px<<", "<<py<<") {"<<i<<"};\n";
	}

	// draw edges
	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(gr); it1 != it2; it1++)
	{
		int s = source(*it1, gr);
		int t = target(*it1, gr);

		assert(v[s] != -1);
		assert(v[t] != -1);

		sprintf(sx, "s%d", s);
		sprintf(sy, "s%d", t);
		assert(s < t);
		
		double bend = -40;
		if(v[s] + 1 == v[t]) bend = 0;
		//else if(v[s] % 2 == 0) bend = -30;

		fout<<"\\draw[line width = 0.02cm, ->, \\colx, bend right = "<< bend <<"] ("<<sx<<") to node {"<< get(get(edge_weight, gr), *it1) <<"} ("<<sy<<");\n";
	}

	draw_footer(fout);

	fout.close();
	return 0;
}
