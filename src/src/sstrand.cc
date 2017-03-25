/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cassert>
#include <cstdio>
#include <map>
#include <iomanip>
#include <fstream>

#include "sstrand.h"

sstrand::sstrand(const string &ch, splice_graph &g)
	: chrm(ch), gr(g)
{
}

sstrand::~sstrand()
{}

int sstrand::build()
{
	build_segments();
	analysis_segments();
	return 0;

	while(remove_inconsistent_strands());
	refine_splice_graph();
	analysis_strand();
	return 0;
}

int sstrand::build_segments()
{
	segments.clear();
	int k = 0;
	int32_t p = -1;
	for(int i = 1; i < gr.num_vertices(); i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		if(vi.lpos != p)
		{
			if(k != 0) segments.push_back(PI(k, i - 1));
			k = i;
		}
		p = vi.rpos;;
	}
	assert(k < gr.num_vertices());
	if(k != 0) segments.push_back(PI(k, gr.num_vertices() - 1));
	return 0;
}

int sstrand::analysis_segments()
{
	for(int i = 0; i < segments.size(); i++)
	{
		analysis_segment(segments[i].first, segments[i].second);
	}
	return 0;
}

int sstrand::analysis_segment(int k1, int k2)
{
	int cnt1 = 0;
	int cnt2 = 0;
	double wrt1 = 0;
	double wrt2 = 0;
	for(int k = k1; k <= k2; k++)
	{
		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.in_edges(k); it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			int t = e->target();
			edge_info ei = gr.get_edge_info(e);
			double w = gr.get_edge_weight(e);
			if(ei.strand == '+') cnt1++;
			if(ei.strand == '-') cnt2++;
			if(ei.strand == '+') wrt1 += w;
			if(ei.strand == '-') wrt2 += w;
		}
		for(tie(it1, it2) = gr.out_edges(k); it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			int t = e->target();
			if(t <= k2) continue;
			edge_info ei = gr.get_edge_info(e);
			double w = gr.get_edge_weight(e);
			if(ei.strand == '+') cnt1++;
			if(ei.strand == '-') cnt2++;
			if(ei.strand == '+') wrt1 += w;
			if(ei.strand == '-') wrt2 += w;
		}
	}

	if(cnt1 == 0 || cnt2 == 0) return 0;

	int32_t p1 = gr.get_vertex_info(k1).lpos;
	int32_t p2 = gr.get_vertex_info(k2).rpos;

	double wrt = wrt1 + wrt2;
	double min = wrt1 < wrt2 ? wrt1 : wrt2;
	double r = min / wrt;

	printf("segment %s:%d-%d cnt %d %d wrt %.0lf %.0lf ratio %.0lf %.0lf %.3lf\n", chrm.c_str(), p1, p2, cnt1, cnt2, wrt1, wrt2, min, wrt, r);

	return 0;
}

bool sstrand::remove_inconsistent_strands()
{
	bool flag = false;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		edge_iterator it1, it2;
		int icnt1 = 0;
		int icnt2 = 0;
		double iwrt1 = 0;
		double iwrt2 = 0;
		VE ive1;
		VE ive2;
		for(tie(it1, it2) = gr.in_edges(i); it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			int t = e->target();
			assert(t == i);
			if(s == 0) continue;
			if(gr.get_vertex_info(s).rpos == gr.get_vertex_info(t).lpos) continue;
			edge_info ei = gr.get_edge_info(e);
			double w = gr.get_edge_weight(e);
			if(ei.strand == '+') icnt1++;
			if(ei.strand == '-') icnt2++;
			if(ei.strand == '+') iwrt1 += w;
			if(ei.strand == '-') iwrt2 += w;
			if(ei.strand == '+') ive1.push_back(e);
			if(ei.strand == '-') ive2.push_back(e);
		}

		int ocnt1 = 0;
		int ocnt2 = 0;
		double owrt1 = 0;
		double owrt2 = 0;
		VE ove1;
		VE ove2;
		for(tie(it1, it2) = gr.out_edges(i); it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			int t = e->target();
			assert(s == i);
			if(t == gr.num_vertices() - 1) continue;
			if(gr.get_vertex_info(s).rpos == gr.get_vertex_info(t).lpos) continue;
			edge_info ei = gr.get_edge_info(e);
			double w = gr.get_edge_weight(e);
			if(ei.strand == '+') ocnt1++;
			if(ei.strand == '-') ocnt2++;
			if(ei.strand == '+') owrt1 += w;
			if(ei.strand == '-') owrt2 += w;
			if(ei.strand == '+') ove1.push_back(e);
			if(ei.strand == '-') ove2.push_back(e);
		}

		if((icnt1 == 0 || icnt2 == 0) && (ocnt1 == 0 || ocnt2 == 0)) continue;

		if(icnt1 >= 1 && iwrt1 < iwrt2 && owrt1 <= owrt2 && ocnt1 == 0)
		{
			// remove edges in ive1
			assert(ive1.size() >= 1);
			flag = true;
			for(int k = 0; k < ive1.size(); k++) gr.remove_edge(ive1[k]);
			printf("remove strands with %lu edges\n", ive1.size());
		}
		else if(icnt2 >= 1 && iwrt2 < iwrt1 && owrt2 <= owrt1 && owrt2 == 0)
		{
			// remove edges in ive2
			assert(ive2.size() >= 1);
			flag = true;
			for(int k = 0; k < ive2.size(); k++) gr.remove_edge(ive2[k]);
			printf("remove strands with %lu edges\n", ive2.size());
		}
		else if(ocnt1 >= 1 && owrt1 < owrt2 && iwrt1 <= iwrt2 && icnt1 == 0)
		{
			// remove edges in ove1
			assert(ove1.size() >= 1);
			flag = true;
			for(int k = 0; k < ove1.size(); k++) gr.remove_edge(ove1[k]);
			printf("remove strands with %lu edges\n", ove1.size());
		}
		else if(ocnt2 >= 1 && owrt2 < owrt1 && iwrt2 <= iwrt1 && iwrt2 == 0)
		{
			// remove edges in ive2
			assert(ove2.size() >= 1);
			flag = true;
			for(int k = 0; k < ove2.size(); k++) gr.remove_edge(ove2[k]);
			printf("remove strands with %lu edges\n", ove2.size());
		}
	}
	return flag;
}


int sstrand::analysis_strand()
{
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		edge_iterator it1, it2;
		int icnt0 = 0;
		int icnt1 = 0;
		int icnt2 = 0;
		double iwrt0 = 0;
		double iwrt1 = 0;
		double iwrt2 = 0;
		for(tie(it1, it2) = gr.in_edges(i); it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			int t = e->target();
			assert(t == i);
			if(s == 0) continue;
			if(gr.get_vertex_info(s).rpos == gr.get_vertex_info(t).lpos) continue;
			edge_info ei = gr.get_edge_info(e);
			double w = gr.get_edge_weight(e);
			if(ei.strand == '.') icnt0++;
			if(ei.strand == '+') icnt1++;
			if(ei.strand == '-') icnt2++;
			if(ei.strand == '.') iwrt0 += w;
			if(ei.strand == '+') iwrt1 += w;
			if(ei.strand == '-') iwrt2 += w;
		}

		int ocnt0 = 0;
		int ocnt1 = 0;
		int ocnt2 = 0;
		double owrt0 = 0;
		double owrt1 = 0;
		double owrt2 = 0;
		for(tie(it1, it2) = gr.out_edges(i); it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			int t = e->target();
			assert(s == i);
			if(t == gr.num_vertices() - 1) continue;
			if(gr.get_vertex_info(s).rpos == gr.get_vertex_info(t).lpos) continue;
			edge_info ei = gr.get_edge_info(e);
			double w = gr.get_edge_weight(e);
			if(ei.strand == '.') ocnt0++;
			if(ei.strand == '+') ocnt1++;
			if(ei.strand == '-') ocnt2++;
			if(ei.strand == '.') owrt0 += w;
			if(ei.strand == '+') owrt1 += w;
			if(ei.strand == '-') owrt2 += w;
		}
		double xi = iwrt1 < iwrt2 ? iwrt1 : iwrt2;
		double xo = owrt1 < owrt2 ? owrt1 : owrt2;
		if(xi <= 0.0 && xo <= 0.1) continue;

		printf("%s:%d-%d icnt = %d %d %d iwrt = %.0lf %.0lf %.0lf %.0lf ocnt = %d %d %d owrt = %.0lf %.0lf %.0lf %.0lf\n", 
				chrm.c_str(), gr.get_vertex_info(i).lpos, gr.get_vertex_info(i).rpos, 
				icnt0, icnt1, icnt2, iwrt0, iwrt1, iwrt2, xi,
				ocnt0, ocnt1, ocnt2, owrt0, owrt1, owrt2, xo);
	}

	return 0;
}

int sstrand::refine_splice_graph()
{
	while(true)
	{
		bool b = false;
		for(int i = 1; i < gr.num_vertices() - 1; i++)
		{
			if(gr.degree(i) == 0) continue;
			if(gr.in_degree(i) >= 1 && gr.out_degree(i) >= 1) continue;
			gr.clear_vertex(i);
			b = true;
		}
		if(b == false) break;
	}
	return 0;
}
