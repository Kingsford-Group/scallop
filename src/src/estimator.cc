#include "estimator.h"

#include <cstdio>
#include <iostream>
#include <cfloat>
#include <algorithm>

estimator::estimator(const splice_graph &g, const vector<path> &p)
	: gr(g), paths(p)
{
}

int estimator::estimate()
{
	if(paths.size() == 0) return 0;

	build_path_indices();
	build_path_lengths();

	print();
	while(true)
	{
		double b = iterate();
		print();
		if(b < 0.01) break;
	}
	printf("\n");

	return 0;
}

int estimator::build_path_indices()
{
	vs.resize(gr.num_vertices());
	for(int i = 0; i < paths.size(); i++)
	{
		for(int k = 0; k < paths[i].v.size(); k++)
		{
			int v = paths[i].v[k];
			assert(v >= 0 && v < gr.num_vertices());
			vs[v].insert(i);
		}
	}
	return 0;
}

int estimator::build_path_lengths()
{
	for(int i = 0; i < paths.size(); i++)
	{
		build_path_length(paths[i]);
	}
	return 0;
}

int estimator::build_path_length(path &p)
{
	p.length = 0;
	for(int i = 0; i < p.v.size(); i++)
	{
		int v = p.v[i];
		if(v == 0) continue;
		if(v == gr.num_vertices() - 1) continue;
		p.length += gr.get_vertex_info(v).length;	
	}
	return 0;
}

double estimator::iterate()
{
	vector<double> weights;
	weights.assign(paths.size(), 0);
	for(int i = 1; i < gr.num_vertices(); i++)
	{
		double ww = gr.get_vertex_weight(i);
		int length = gr.get_vertex_info(i).length;
		set<int> &s = vs[i];
		double sum = 0;
		double pseudo = 0.01;
		for(set<int>::iterator it = s.begin(); it != s.end(); it++)
		{
			int p = (*it);
			sum += paths[p].abd + pseudo;
		}
		for(set<int>::iterator it = s.begin(); it != s.end(); it++)
		{
			int p = (*it);
			double ratio = (paths[p].abd + pseudo) / sum;
			weights[p] += ww * length * ratio;
		}
	}

	double error = 0;
	for(int i = 0; i < paths.size(); i++)
	{
		double abd = weights[i] / paths[i].length;
		error += fabs(abd - paths[i].abd) / paths[i].abd;
		paths[i].abd = abd;
	}
	return error / paths.size();
}

int estimator::print()
{
	printf("abundance: ");
	for(int i = 0; i < paths.size(); i++) printf("%d:%.2lf ", i, paths[i].abd);
	printf("\n");
	return 0;
}
