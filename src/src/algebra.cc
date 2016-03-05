#include "algebra.h"

#include <iomanip>
#include <cstdio>
#include <cfloat>

algebra::algebra()
{}

algebra::~algebra()
{}

int algebra::solve(int n, int m)
{
	srand(time(0));
	simulate(n, m);
	write_splice_graph(gr, "sgraph.gr");
	draw_splice_graph(gr, "sgraph.tex");
	process();
	return 0;
}

int algebra::solve(const string &file)
{
	build_splice_graph(gr, file);
	process();
	return 0;
}

int algebra::process()
{
	get_edge_indices(gr, i2e, e2i);
	build_path_edge_matrix();
	pem0 = pem;
	sample_valid_paths();

	print();
	gaussian_elimination();
	print();

	check();
	return 0;
}

int algebra::simulate(int n, int m)
{
	while(true)
	{
		simulate_splice_graph(gr, n, m);
		bool b = check_fully_connected(gr);
		if(b == true) break;
	}
	return 0;
}

int algebra::sample_valid_paths()
{
	for(int i = 0; i < 100; i++)
	{
		sample_paths();
		bool b = check_edges_cover();
		if(b == true) return 0;
	}
	return 0;
}

int algebra::sample_paths()
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

int algebra::check_edges_cover()
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

int algebra::build_path_edge_matrix()
{
	vector<int> v;
	build_path_edge_matrix(0, v);
	return 0;
}

int algebra::build_path_edge_matrix(int s, vector<int> &v)
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

int algebra::gaussian_elimination()
{
	rank = pem.size();
	for(int i = 0; i < pem.size(); i++)
	{
		int f = gaussian_elimination(i);
		if(f == -1) 
		{
			rank = i;
			return 0;
		}
	}
	return 0;
}

int algebra::gaussian_elimination(int r)
{
	int f = choose_main_element(r);
	if(f == -1) return -1;
	normalize_row(r);
	eliminate_column(r);
	return 0;
}

int algebra::eliminate_column(int r)
{
	for(int k = 0; k < pem.size(); k++)
	{
		if(k == r) continue;
		if(fabs(pem[k][r]) < SMIN) continue;
		add_to_row(r, k, 0 - pem[k][r]);
	}
	return 0;
}

int algebra::add_to_row(int s, int t, double c)
{
	for(int i = 0; i < pem[0].size(); i++)
	{
		pem[t][i] += pem[s][i] * c;
	}
	return 0;
}

int algebra::normalize_row(int r)
{
	assert(fabs(pem[r][r]) >= SMIN);
	for(int i = 0; i < r; i++) assert(fabs(pem[r][i]) < SMIN);

	for(int i = r + 1; i < pem[0].size(); i++)
	{
		pem[r][i] = pem[r][i] / pem[r][r];
	}
	pem[r][r] = 1;
	return 0;
}

int algebra::choose_main_row(int r)
{
	if(pem[r][r] != 0) return 0;
	int k = -1;
	for(int i = r + 1; i < pem.size(); i++)
	{
		if(pem[i][r] != 0)
		{
			k = i;
			break;
		}
	}
	if(k == -1) return -1;
	vector<double> v = pem[r];
	pem[r] = pem[k];
	pem[k] = v;
	return 0;
}

int algebra::choose_main_column(int r)
{
	if(pem[r][r] != 0) return 0;
	int k = -1;
	for(int i = r + 1; i < pem[0].size(); i++)
	{
		if(pem[r][i] != 0)
		{
			k = i;
			break;
		}
	}
	if(k == -1) return -1;
	vector<double> v;
	for(int i = 0; i < pem.size(); i++) v.push_back(pem[i][r]);
	for(int i = 0; i < pem.size(); i++) pem[i][r] = pem[i][k];
	for(int i = 0; i < pem.size(); i++) pem[i][k] = v[i];
	return 0;
}

int algebra::choose_main_element(int r)
{
	int ii = -1;
	int jj = -1;
	for(int j = r; j < pem[0].size(); j++)
	{
		for(int i = r; i < pem.size(); i++)
		{
			if(fabs(pem[i][j]) > SMIN)
			{
				ii = i;
				jj = j;
				break;
			}
		}
	}
	if(ii == -1) return -1;
	if(ii != r) exchange_row(r, ii);
	if(jj != r) exchange_column(r, jj);
	return 0;
}

int algebra::exchange_row(int s, int t)
{
	vector<double> v = pem[s];
	pem[s] = pem[t];
	pem[t] = v;
	return 0;
}

int algebra::exchange_column(int s, int t)
{
	vector<double> v;
	for(int i = 0; i < pem.size(); i++) v.push_back(pem[i][s]);
	for(int i = 0; i < pem.size(); i++) pem[i][s] = pem[i][t];
	for(int i = 0; i < pem.size(); i++) pem[i][t] = v[i];
	return 0;
}

int algebra::check()
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


int algebra::print()
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
