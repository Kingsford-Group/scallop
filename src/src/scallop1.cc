#include "scallop1.h"
#include <cstdio>

scallop1::scallop1(splice_graph &gr)
	: assembler(gr)
{}

scallop1::~scallop1()
{}

int scallop1::assemble()
{
	smooth_weights();
	iterate();
	return 0;
}

int scallop1::iterate()
{
	while(true)
	{
		path p0 = compute_maximum_forward_path();
		
		path max_qx, max_qy;
		double max_gain = 0;
		int max_index = -1;
		for(int k = 0; k < paths.size(); k++)
		{
			path &px = paths[k];
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
			paths.push_back(p0);
		}
		else
		{
			printf("MAX GAIN = %.4lf\n", max_gain);
			increase_path(paths[max_index]);
			decrease_path(max_qx);
			decrease_path(max_qy);
			paths[max_index] = max_qx;
			paths.push_back(max_qy);
		}

		if(max_gain + p0.abd < 0.1) break;
	}
	return 0;
}

int scallop1::resolve(const path &px, const path &py, path &qx, path &qy) const
{
	assert(px.abd >= py.abd);
	vector< vector<int> > vv;
	vv.resize(2);
	vv[0] = px.index(gr.num_vertices());
	vv[1] = py.index(gr.num_vertices());

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
		if(t == gr.num_vertices() - 1) break;
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
		if(t == gr.num_vertices() - 1) break;
		s = t;
	}

	qx.abd = py.abd;
	qy.abd = py.abd;
	return 0;
}

