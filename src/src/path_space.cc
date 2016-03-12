#include "path_space.h"

#include <iomanip>
#include <cstdio>
#include <cfloat>
#include <cstdio>
#include <cmath>

path_space::path_space()
{
	rank = -1;
}

path_space::~path_space()
{}

int path_space::solve(int n, int m)
{
	srand(time(0));
	simulate(n, m);
	write_splice_graph(gr, "sgraph.gr");
	draw_splice_graph(gr, "sgraph.tex");
	process();
	return 0;
}

int path_space::solve(const string &file)
{
	build_splice_graph(gr, file);
	process();
	return 0;
}

int path_space::process()
{
	get_edge_indices(gr, i2e, e2i);
	build_path_edge_matrix();
	pem0 = pem;
	sample_valid_paths();

	algebra ab(pem);
	rank = ab.full_eliminate();
	pem = ab.mat;
	print();

	check();
	return 0;
}

int path_space::simulate(int n, int m)
{
	while(true)
	{
		simulate_splice_graph(gr, n, m);
		bool b = check_fully_connected(gr);
		if(b == true) break;
	}
	return 0;
}

int path_space::sample_valid_paths()
{
	for(int i = 0; i < 100; i++)
	{
		sample_paths();
		bool b = check_edges_cover();
		if(b == true) return 0;
	}
	return 0;
}

int path_space::sample_paths()
{
	vector<int> v;
	for(int i = 0; i < pem0.size(); i++) v.push_back(i);

	pem.clear();
	int r = num_edges(gr) - num_vertices(gr) + 1;
	for(int i = 0; i < r; i++)
	{
		int p = rand() % (v.size() - i);
		pem.push_back(pem0[v[p]]);
		int x = v[p];
		v[p] = v[v.size() - i - 1];
		v[v.size() - i - 1] = x;
	}

	return 0;
}

int path_space::check_edges_cover()
{
	for(int j = 0; j < pem[0].size(); j++)
	{
		bool b = false;
		for(int i = 0; i < pem.size(); i++)
		{
			if(fabs(pem[i][j]) > SMIN) b = true;
		}
		if(b == false) return false;
	}
	return true;
}

int path_space::build_path_edge_matrix()
{
	vector<int> v;
	build_path_edge_matrix(0, v);
	return 0;
}

int path_space::build_path_edge_matrix(int s, vector<int> &v)
{
	if(s == num_vertices(gr) - 1)
	{
		vector<double> p;
		p.resize(num_edges(gr));
		for(int i = 0; i < v.size(); i++)
		{
			p[v[i]] = 1;
		}
		pem.push_back(p);
		return 0;
	}

	out_edge_iterator it1, it2;
	for(tie(it1, it2) = out_edges(s, gr); it1 != it2; it1++)
	{
		int t = target(*it1, gr);
		v.push_back(e2i[*it1]);
		build_path_edge_matrix(t, v);
		v.pop_back();
	}
	return 0;
}

int path_space::check()
{
	for(int j = rank; j < pem[0].size(); j++)
	{
		for(int i = 0; i < pem.size(); i++)
		{
			bool b1 = fabs(pem[i][j]) < SMIN ? true : false;
			bool b2 = fabs(fabs(pem[i][j]) - 1) < SMIN ? true : false;
			if(b1 || b2) continue;
			printf("EXAMPLE\n");
		}
	}
	return 0;
}	


int path_space::print()
{
	printf("     ");
	for(int j = 0; j < pem[0].size(); j++)
	{
		printf("[%2d] ", j + 1);
	}
	printf("\n");

	for(int i = 0; i < pem.size(); i++)
	{
		printf("[%2d] ", i + 1);
		for(int j = 0; j < pem[i].size(); j++)
		{
			if(fabs(pem[i][j]) < SMIN) printf("%4.1lf ", 0.0);
			else printf("%4.1lf ", pem[i][j]);
		}
		printf("\n");
	}
	int m = num_edges(gr);
	int n = num_vertices(gr);
	int r = m - n + 2;
	printf("edges = %d, vertices = %d, rank = %d\n\n", m, n, rank);
	return 0;
}
