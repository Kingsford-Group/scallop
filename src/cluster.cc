/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "cluster.h"
#include "config.h"
#include <cassert>
#include <algorithm>

cluster::cluster(const vector<transcript> &v)
	:trs(v)
{
}

int cluster::solve()
{
	split_with_num_exons();
	build_graph();
	clustering();
	return 0;
}

int cluster::split_with_num_exons()
{
	for(int i = 0; i < trs.size(); i++)
	{
		int e = trs[i].exons.size();
		if(trs[i].strand == '-') e = 0 - e;
		if(miv.find(e) == miv.end())
		{
			vector<int> v;
			v.push_back(i);
			miv.insert(PIV(e, v));
		}
		else
		{
			miv[e].push_back(i);
		}
	}
	return 0;
}

int cluster::build_graph()
{
	gr.clear();
	for(int i = 0; i < trs.size(); i++) gr.add_vertex();
	
	for(MIV::iterator it = miv.begin(); it != miv.end(); it++)
	{
		int e = it->first;
		vector<int> &v = it->second;
		for(int i = 0; i < v.size(); i++)
		{
			for(int j = 0; j < i; j++)
			{
				if(verify_equal(v[i], v[j])) gr.add_edge(v[i], v[j]);
				//if(verify_subset(v[i], v[j])) gr.add_edge(v[i], v[j]);
				//if(verify_subset(v[j], v[i])) gr.add_edge(v[i], v[j]);
			}
		}
	}
	return 0;
}

int cluster::clustering()
{
	cct.clear();
	vector< set<int> > cc = gr.compute_connected_components();
	for(int i = 0; i < cc.size(); i++)
	{
		set<int> &s = cc[i];

		int k = -1;				// transcript index with max-capacity (or longest)
		int32_t lpos = -1;		// leftmost position
		int32_t rpos = -1;		// rightmost position
		double cov = 0;			// sum of coverage
		vector<PI32> chain;		// intron chain
		vector<int> cnts;		// counts for each intron
		for(set<int>::iterator it = s.begin(); it != s.end(); it++)
		{
			transcript &t = trs[*it];
			cov += t.coverage;
			PI32 p = t.get_bounds();
			if(lpos == -1 || p.first < lpos) lpos = p.first;
			if(rpos == -1 || p.second > rpos) rpos = p.second;
			if(k == -1 || t.coverage > trs[k].coverage) k = (*it);
		}

		if(k == -1) continue;
		transcript cc = trs[k];
		cc.coverage = cov;
		cc.exons[0].first = lpos;
		cc.exons[cc.exons.size() - 1].second = rpos;

		cct.push_back(cc);

		if(verbose >= 2) printf("cct: %lu transcripts, exons = %lu, max-coverage = %.2lf, sum-coverage = %.2lf\n", 
				s.size(), trs[k].exons.size(), trs[k].coverage, cov);
	}

	return 0;
}

bool cluster::verify_equal(int x, int y)
{
	transcript &xx = trs[x];
	transcript &yy = trs[y];
	if(xx.strand != yy.strand) return false;
	if(xx.exons.size() != yy.exons.size()) return false;

	PI32 px = xx.get_bounds();
	PI32 py = yy.get_bounds();
	if(xx.exons.size() == 1)
	{
		if(px.first <= py.first)
		{
			if(px.second >= py.second) return true;
			if(px.second <= py.first) return false;
			double r1 = (px.second - py.first) * 1.0 / (px.second - px.first);
			double r2 = (px.second - py.first) * 1.0 / (py.second - py.first);
			double r = r1 > r2 ? r1 : r2;
			if(r >= min_cluster_single_exon_ratio) return true;
			else return false;
		}
		else
		{
			if(py.second >= px.second) return true;
			if(py.second <= px.first) return false;
			double r1 = (py.second - px.first) * 1.0 / (py.second - py.first);
			double r2 = (py.second - px.first) * 1.0 / (px.second - px.first);
			double r = r1 > r2 ? r1 : r2;
			if(r >= min_cluster_single_exon_ratio) return true;
			else return false;
		}
	}

	if(fabs(px.first - py.first) > max_cluster_boundary_distance) return false;
	if(fabs(px.second - py.second) > max_cluster_boundary_distance) return false;

	vector<PI32> ix = xx.get_intron_chain();
	vector<PI32> iy = yy.get_intron_chain();
	assert(ix.size() == iy.size());
	for(int i = 0; i < ix.size(); i++)
	{
		double f1 = fabs(ix[i].first - iy[i].first);
		double f2 = fabs(ix[i].second - iy[i].second);
		if(f1 + f2 > max_cluster_intron_distance) return false;
	}
	return true;
}

bool cluster::verify_subset(int x, int y)
{
	transcript &xx = trs[x];
	transcript &yy = trs[y];
	if(xx.strand != yy.strand) return false;
	if(xx.exons.size() >= yy.exons.size()) return false;
	if(xx.exons.size() == 1) return false;
	if(xx.coverage >= yy.coverage + 0.000001) return false;

	PI32 px = xx.get_bounds();
	PI32 py = yy.get_bounds();
	vector<PI32> ix = xx.get_intron_chain();
	vector<PI32> iy = yy.get_intron_chain();

	if(px.first < py.first) return false;
	if(px.second > py.second) return false;

	int n = xx.exons.size();
	for(int k = 0; k <= yy.exons.size() - xx.exons.size(); k++)
	{
		if(xx.exons[0].first < yy.exons[k].first) continue;
		if(xx.exons[n - 1].second > yy.exons[k + n - 1].second) continue;

		bool b = true;
		for(int i = 0; i < ix.size(); i++)
		{
			double f1 = fabs(ix[i].first - iy[i + k].first);
			double f2 = fabs(ix[i].second - iy[i + k].second);
			if(f1 + f2 > max_cluster_intron_distance) b = false;
			if(b == false) break;
		}
		if(b == false) continue;

		return true;
	}
}

int cluster::print()
{
	for(MIV::iterator it = miv.begin(); it != miv.end(); it++)
	{
		printf("there are %lu transcripts with %d exons\n", it->second.size(), it->first);
	}
	return 0;
}
