#include "assembler.h"
#include "draw.h"
#include "lpsolver.h"
#include "fibonacci_heap.h"

#include <iomanip>
#include <cfloat>

assembler::assembler(splice_graph &g)
	: gr(g)
{}

assembler::~assembler()
{}

int assembler::solve()
{
	update_weights();
	greedy();
	iterate();
	print();
	return 0;
}

int assembler::update_weights()
{
	lpsolver lp(gr);
	lp.solve();
	return 0;
}

int assembler::greedy()
{
	MED med;
	backup_edge_weights(med);
	while(true)
	{
		path p = compute_maximum_forward_path();
		decrease_path(p);
		if(p.v.size() <= 1) break;
		paths0.push_back(p);
		if(p.abd < 1) break;
	}
	recover_edge_weights(med);
	return 0;
}

int assembler::iterate()
{
	MED med;
	backup_edge_weights(med);
	while(true)
	{
		path p0 = compute_maximum_forward_path();
		
		path max_qx, max_qy;
		double max_gain = 0;
		int max_index = -1;
		for(int k = 0; k < paths1.size(); k++)
		{
			path &px = paths1[k];
			if(px.abd <= p0.abd) continue;
			add_backward_path(px);
			path py = compute_maximum_path();
			remove_backward_path(px);

			//if(2 * py.abd - p0.abd - px.abd < 0) continue; 

			path qx, qy;
			resolve(px, py, qx, qy);
			assert(qx.abd == py.abd);
			assert(qy.abd == py.abd);

			increase_path(px);

			decrease_path(qx);
			double ay = compute_bottleneck_weight(qy);
			increase_path(qx);

			decrease_path(qy);
			double bx = compute_bottleneck_weight(qx);
			increase_path(qy);

			decrease_path(px);

			double a = qx.abd + ay - px.abd - p0.abd;
			double b = bx + qy.abd - px.abd - p0.abd;
			double w = (a > b) ? a : b;

			if(w <= max_gain + 0.01) continue;

			assert(qx.abd == qy.abd);
			
			if(a > b) qy.abd = ay;
			else qx.abd = bx;

			max_index = k;
			max_gain = w;
			max_qx = qx;
			max_qy = qy;
		}

		if(max_index == -1)
		{
			if(p0.v.size() < 2) break;
			decrease_path(p0);
			paths1.push_back(p0);
		}
		else
		{
			printf("MAX GAIN = %.4lf\n", max_gain);
			increase_path(paths1[max_index]);
			decrease_path(max_qx);
			decrease_path(max_qy);
			paths1[max_index] = max_qx;
			paths1.push_back(max_qy);
		}

		if(max_gain + p0.abd < 1) break;
	}
	recover_edge_weights(med);
	return 0;
}

path assembler::compute_maximum_forward_path() const
{
	path p;
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
			//if(s >= i) continue;
			double xw = get(get(edge_weight, gr), *it1);
			double ww = xw < table[s] ? xw : table[s];
			if(ww >= max_abd)
			{
				max_abd = ww;
				max_idx = s;
			}
		}

		back[i] = max_idx;
		table[i] = max_abd < abd ? max_abd : abd;
	}

	if(table[n - 1] <= 0) return p;

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

	return p;
}

path assembler::compute_maximum_path() const
{
	path p;
	fibonacci_heap f;						// fibonacci heap
	vector<handle_t> handles;		// handles for nodes
	vector<bool> reached;			// whether each node is reached
	vector<bool> visited;			// whether each node is visited
	vector<int> backptr;			// backward pointers

	handles.resize(num_vertices(gr));
	reached.resize(num_vertices(gr), false);
	visited.resize(num_vertices(gr), false);
	backptr.resize(num_vertices(gr), -1);

	handle_t h = f.push(fnode(0, DBL_MAX));
	handles[0] = h;
	reached[0] = true;
	backptr[0] = -1;

	while(f.empty() == false)
	{
		fnode fx = f.top();
		visited[fx.v] = true;

		if(fx.v == num_vertices(gr) - 1) break;
		
		f.pop();
		out_edge_iterator it1, it2;
		for(tie(it1, it2) = out_edges(fx.v, gr); it1 != it2; it1++)
		{
			int y = target(*it1, gr);
			if(visited[y] == true) continue;

			double xw = get(get(edge_weight, gr), *it1);
			double ww = xw < fx.w ? xw : fx.w;
			if(reached[y] == false) 
			{
				h = f.push(fnode(y, ww));
				handles[y] = h;
				reached[y] = true;
				backptr[y] = fx.v;
			}
			else
			{
				double yw = (*handles[y]).w;
				if(ww > yw) 
				{
					f.increase(handles[y], fnode(y, ww));
					backptr[y] = fx.v;
				}
			}
		}
	}

	if(f.empty() == true) return p;

	p.abd = f.top().w;
	p.v.clear();
	int b = num_vertices(gr) - 1;
	while(true)
	{
		p.v.push_back(b);
		if(b == 0) break;
		b = backptr[b];
		assert(b != -1);
	}
	reverse(p.v.begin(), p.v.end());

	return p;
}

int assembler::decrease_path(const path &p)
{
	if(p.v.size() < 2) return 0;
	for(int i = 0; i < p.v.size() - 1; i++)
	{
		PEB e = edge(p.v[i], p.v[i + 1], gr);
		assert(e.second == true);
		double w0 = get(get(edge_weight, gr), e.first);
		double w1 = w0 - p.abd;
		
		assert(w1 >= -0.000001);
		if(w1 <= 0) w1 = 0;
		put(get(edge_weight, gr), e.first, w1);
	}
	return 0;
}

int assembler::increase_path(const path &p)
{
	if(p.v.size() < 2) return 0;
	for(int i = 0; i < p.v.size() - 1; i++)
	{
		PEB e = edge(p.v[i], p.v[i + 1], gr);
		assert(e.second == true);
		double w0 = get(get(edge_weight, gr), e.first);
		double w1 = w0 + p.abd;
		put(get(edge_weight, gr), e.first, w1);
	}
	return 0;
}


int assembler::add_backward_path(const path &p)
{
	if(p.v.size() < 2) return 0;
	for(int i = 0; i < p.v.size() - 1; i++)
	{
		PEB e = add_edge(p.v[i + 1], p.v[i], gr);
		assert(e.second == true);
		put(get(edge_weight, gr), e.first, p.abd);
	}
	return 0;
}

int assembler::remove_backward_path(const path &p)
{
	if(p.v.size() < 2) return 0;
	for(int i = 0; i < p.v.size() - 1; i++)
	{
		remove_edge(p.v[i + 1], p.v[i], gr);
	}
	return 0;
}

double assembler::compute_bottleneck_weight(const path &p) const
{
	double ww = DBL_MAX;
	for(int i = 0; i < p.v.size() - 1; i++)
	{
		PEB e = edge(p.v[i], p.v[i + 1], gr);
		assert(e.second == true);
		double w = get(get(edge_weight, gr), e.first);
		if(w < ww) ww = w;
	}
	return ww;
}

int assembler::resolve(const path &px, const path &py, path &qx, path &qy) const
{
	assert(px.abd >= py.abd);
	vector< vector<int> > vv;
	vv.resize(2);
	vv[0] = px.index(num_vertices(gr));
	vv[1] = py.index(num_vertices(gr));

	int k = 0;
	int s = 0;
	qx.clear();
	qx.v.push_back(s);
	while(true)
	{
		int t = vv[k][s];
		if(vv[1 - k][t] == s)
		{
			k = 1 - k;
			t = vv[k][s];
		}
		qx.v.push_back(t);
		if(t == num_vertices(gr) - 1) break;
		s = t;
	}

	k = 1;
	s = 0;
	qy.clear();
	qy.v.push_back(s);
	while(true)
	{
		int t = vv[k][s];
		if(vv[1 - k][t] == s)
		{
			k = 1 - k;
			t = vv[k][s];
		}
		qy.v.push_back(t);
		if(t == num_vertices(gr) - 1) break;
		s = t;
	}

	qx.abd = py.abd;
	qy.abd = py.abd;
	return 0;
}

int assembler::backup_edge_weights(MED &med) const
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

int assembler::recover_edge_weights(const MED &med)
{
	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(gr); it1 != it2; it1++)
	{
		MED::const_iterator it = med.find(*it1);
		put(get(edge_weight, gr), *it1, it->second);
	}
	return 0;
}

int assembler::print() const
{
	// print paths0
	double w0 = 0.0;
	for(int i = 0; i < paths0.size(); i++)
	{
		printf("greedy  ");
		paths0[i].print(i);
		w0 += paths0[i].abd;	
	}

	// print paths1
	double w1 = 0.0;
	for(int i = 0; i < paths1.size(); i++)
	{
		printf("iterate ");
		paths1[i].print(i);
		w1 += paths1[i].abd;	
	}
	printf("greedy  summary: %ld paths, %.2lf total abundance\n", paths0.size(), w0);
	printf("iterate summary: %ld paths, %.2lf total abundance\n", paths1.size(), w1);

	return 0;
}

int assembler::draw(const string &file) const
{
	ofstream fout(file.c_str());
	if(fout.fail())
	{
		printf("open file %s error.\n", file.c_str());
		return 0;
	}

	draw_header(fout);

	double len = 2.4;
	fout<<"\\def\\len{"<<len<<"cm}\n";

	// draw vertices
	char sx[1024];
	char sy[1024];
	for(int i = 0; i < num_vertices(gr); i++)
	{
		sprintf(sx, "s%d", i);
		fout.precision(1);
		fout<<fixed;
		fout<<"\\node[mycircle, \\colx, draw, label = below:{";
		fout<< get(get(vertex_weight, gr), i) << "," << get(get(vertex_stddev, gr), i);
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
		fout<< get(get(edge_weight, gr), *it1) <<",";
		fout<< get(get(edge_stddev, gr), *it1) <<"} ("<<sy<<");\n";
	}

	draw_footer(fout);

	fout.close();
	return 0;
}

