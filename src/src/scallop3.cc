#include "scallop3.h"
#include <cstdio>

scallop3::scallop3(splice_graph &gr)
	: assembler(gr)
{}

scallop3::~scallop3()
{}

int scallop3::assemble()
{
	smooth_weights();
	init_super_edges();
	reconstruct_splice_graph();
	return 0;
}

int scallop3::init_super_edges()
{
	mev.clear();
	edge_iterator it1, it2;
	vector<int> v;
	for(tie(it1, it2) = edges(gr); it1 != it2; it1++)
	{
		mev.insert(PEV(*it1, v));
	}
	return 0;
}

int scallop3::reconstruct_splice_graph()
{
	while(true)
	{
		bool flag = false;
		for(int i = 0; i < num_vertices(gr); i++)
		{
			bool b = decompose_trivial_vertex(i);
			if(b == true) flag = true;
		}
		if(flag == false) break;
	}
	return 0;
}

bool scallop3::decompose_trivial_vertex(int x)
{
	int id = in_degree(x, gr);
	int od = out_degree(x, gr);

	if(id <= 0 || od <= 0) return false;
	if(id >= 2 && od >= 2) return false;
	//if(id <= 1 && od <= 1) return false;

	in_edge_iterator it1, it2;
	out_edge_iterator ot1, ot2;
	for(tie(it1, it2) = in_edges(x, gr); it1 != it2; it1++)
	{

		for(tie(ot1, ot2) = out_edges(x, gr); ot1 != ot2; ot1++)
		{
			int s = source(*it1, gr);
			int t = target(*ot1, gr);

			double w1 = get(get(edge_weight, gr), *it1);
			double a1 = get(get(edge_stddev, gr), *it1);
			double w2 = get(get(edge_weight, gr), *ot1);
			double a2 = get(get(edge_stddev, gr), *ot1);

			double w = w1 < w2 ? w1 : w2;
			double a = w1 < w2 ? a1 : a2;

			PEB p = add_edge(s, t, gr);
			put(get(edge_weight, gr), p.first, w);
			put(get(edge_stddev, gr), p.first, a);

			assert(mev.find(*it1) != mev.end());
			assert(mev.find(*ot1) != mev.end());

			vector<int> v = mev[*it1];
			v.push_back(x);
			v.insert(v.end(), mev[*ot1].begin(), mev[*ot1].end());

			mev.insert(PEV(p.first, v));
		}
	}
	clear_vertex(x, gr);
	return true;
}

int scallop3::print()
{
	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(gr); it1 != it2; it1++)
	{
		int s = source(*it1, gr);
		int t = target(*it1, gr);
		double w = get(get(edge_weight, gr), *it1);
		printf("edge: (%5d -> %5d), weight = %.0lf\n", s, t, w);
	}
	return 0;
}

bool scallop3::decide_nested()
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
				if(check_directed_path(gr, i, k) == false) continue;
				if(check_directed_path(gr, k, j) == false) continue;
				out_edge_iterator it3, it4;
				for(tie(it3, it4) = out_edges(k, gr); it3 != it4; it3++)
				{
					int l = target(*it3, gr);
					assert(l > k);
					if(l <= j) continue;
					
					if(check_directed_path(gr, j, l) == false) continue;

					return false;
				}
			}
		}
	}
	return true;
}


