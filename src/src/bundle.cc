#include <cassert>
#include <cstdio>
#include <map>
#include <iomanip>
#include <fstream>

#include "bundle.h"
#include "binomial.h"
#include "region.h"
#include "config.h"
#include "util.h"
#include "undirected_graph.h"

bundle::bundle(const bundle_base &bb)
	: bundle_base(bb)
{
}

bundle::~bundle()
{}

int bundle::build()
{
	//compute_strand();
	check_left_ascending();

	build_junctions();

	bool bb = identify_extra_boundary;

	if(identify_extra_boundary == true)
	{
		build_regions(5);
		build_partial_exons();

		identify_extra_boundary = false;

		build_regions(0);
		build_partial_exons();

		identify_extra_boundary = true;
	}
	else
	{
		build_regions(0);
		build_partial_exons();
	}

	build_partial_exon_map();
	link_partial_exons();
	build_splice_graph();

	// revise splice graph
	//join_single_exon_transcripts();

	remove_small_edges();
	refine_splice_graph();

	remove_inner_vertices();
	remove_inner_start_boundaries();
	remove_inner_end_boundaries();

	if(extend_isolated_boundary == true)
	{
		extend_isolated_start_boundaries();
		extend_isolated_end_boundaries();
	}

	build_hyper_edges2();
	//gr.draw("gr.tex");

	return 0;
}

int bundle::compute_strand()
{
	int n0 = 0, np = 0, nq = 0;
	for(int i = 0; i < hits.size(); i++)
	{
		if(hits[i].xs == '.') n0++;
		if(hits[i].xs == '+') np++;
		if(hits[i].xs == '-') nq++;
	}

	if(np > nq) strand = '+';
	if(np < nq) strand = '-';

	return 0;
}

int bundle::check_left_ascending()
{
	for(int i = 1; i < hits.size(); i++)
	{
		int32_t p1 = hits[i - 1].pos;
		int32_t p2 = hits[i].pos;
		assert(p1 <= p2);
	}
	return 0;
}

int bundle::check_right_ascending()
{
	for(int i = 1; i < hits.size(); i++)
	{
		int32_t p1 = hits[i - 1].rpos;
		int32_t p2 = hits[i].rpos;
		assert(p1 <= p2);
	}
	return 0;
}

int bundle::build_junctions()
{
	map<int64_t, int> m;
	vector<int64_t> v;
	for(int i = 0; i < hits.size(); i++)
	{
		hits[i].get_splice_positions(v);
		if(v.size() == 0) continue;

		//hits[i].print();
		for(int k = 0; k < v.size(); k++)
		{
			int64_t p = v[k];
			//printf(" %d-%d\n", low32(p), high32(p));
			if(m.find(p) == m.end()) m.insert(pair<int64_t, int>(p, 1));
			else m[p]++;
		}
	}

	map<int64_t, int>::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		if(it->second < min_splice_boundary_hits) continue;
		junctions.push_back(junction(it->first, it->second));
	}
	return 0;
}

int bundle::build_regions(int count)
{
	MPI s;
	s.insert(PI(lpos, START_BOUNDARY));
	s.insert(PI(rpos, END_BOUNDARY));
	for(int i = 0; i < junctions.size(); i++)
	{
		junction &jc = junctions[i];

		double ave, dev;
		evaluate_rectangle(mmap, jc.lpos, jc.rpos, ave, dev);
		if(jc.count < count && ave >= count) continue;

		int32_t l = jc.lpos;
		int32_t r = jc.rpos;

		if(s.find(l) == s.end()) s.insert(PI(l, LEFT_SPLICE));
		else if(s[l] == RIGHT_SPLICE) s[l] = LEFT_RIGHT_SPLICE;

		if(s.find(r) == s.end()) s.insert(PI(r, RIGHT_SPLICE));
		else if(s[r] == LEFT_SPLICE) s[r] = LEFT_RIGHT_SPLICE;
	}

	for(int i = 0; i < pexons.size(); i++)
	{
		partial_exon &p = pexons[i];
		if(s.find(p.lpos) != s.end()) s.insert(PI(p.lpos, p.ltype));
		if(s.find(p.rpos) != s.end()) s.insert(PI(p.rpos, p.rtype));
	}

	vector<PPI> v(s.begin(), s.end());
	sort(v.begin(), v.end());

	regions.clear();
	for(int k = 0; k < v.size() - 1; k++)
	{
		int32_t l = v[k].first;
		int32_t r = v[k + 1].first;
		int ltype = v[k].second; 
		int rtype = v[k + 1].second; 

		if(ltype == LEFT_RIGHT_SPLICE) ltype = RIGHT_SPLICE;
		if(rtype == LEFT_RIGHT_SPLICE) rtype = LEFT_SPLICE;

		regions.push_back(region(l, r, ltype, rtype, &mmap, &imap));
	}

	return 0;
}

int bundle::build_partial_exons()
{
	pexons.clear();
	for(int i = 0; i < regions.size(); i++)
	{
		region &r = regions[i];
		pexons.insert(pexons.end(), r.pexons.begin(), r.pexons.end());
	}
	return 0;
}

int bundle::build_partial_exon_map()
{
	pmap.clear();
	for(int i = 0; i < pexons.size(); i++)
	{
		partial_exon &p = pexons[i];
		pmap += make_pair(ROI(p.lpos, p.rpos), i + 1);
	}
	return 0;
}

int bundle::locate_left_partial_exon(int32_t x)
{
	SIMI it = pmap.find(ROI(x, x + 1));
	if(it == pmap.end()) return -1;
	assert(it->second >= 1);
	assert(it->second <= pexons.size());

	int k = it->second - 1;
	int32_t p1 = lower(it->first);
	int32_t p2 = upper(it->first);
	assert(p2 >= x);
	assert(p1 <= x);

	if(x - p1 > min_flank_length && p2 - x < min_flank_length) k++;
	if(k >= pexons.size()) return -1;
	return k;
}

int bundle::locate_right_partial_exon(int32_t x)
{
	SIMI it = pmap.find(ROI(x - 1, x));
	if(it == pmap.end()) return -1;
	assert(it->second >= 1);
	assert(it->second <= pexons.size());

	int k = it->second - 1;
	int32_t p1 = lower(it->first);
	int32_t p2 = upper(it->first);
	assert(p1 < x);
	assert(p2 >= x);

	if(p2 - x > min_flank_length && x - p1 <= min_flank_length) k--;
	return k;
}

int bundle::build_hyper_edges1()
{
	hs.clear();
	for(int i = 0; i < hits.size(); i++)
	{
		hit &h = hits[i];
		if((h.flag & 0x4) >= 1) continue;

		vector<int64_t> v;
		h.get_matched_intervals(v);
		if(v.size() == 0) continue;

		set<int> sp;
		for(int k = 0; k < v.size(); k++)
		{
			int32_t p1 = high32(v[k]);
			int32_t p2 = low32(v[k]);

			int k1 = locate_left_partial_exon(p1);
			int k2 = locate_right_partial_exon(p2);
			if(k1 < 0 || k2 < 0) continue;

			for(int j = k1; j <= k2; j++) sp.insert(j);
		}

		if(sp.size() <= 1) continue;
		hs.add_node_list(sp);
	}

	return 0;
}

int bundle::build_hyper_edges2()
{
	//sort(hits.begin(), hits.end(), hit_compare_by_name);
	sort(hits.begin(), hits.end());

	/*
	printf("----------------------\n");
	for(int k = 9; k < hits.size(); k++) hits[k].print();
	printf("======================\n");
	*/

	hs.clear();

	string qname;
	int hi = -2;
	vector<int> sp1;
	for(int i = 0; i < hits.size(); i++)
	{
		hit &h = hits[i];
		
		/*
		printf("sp1 = ( ");
		printv(sp1);
		printf(")\n");
		h.print();
		*/

		if(h.qname != qname || h.hi != hi)
		{
			set<int> s(sp1.begin(), sp1.end());
			if(s.size() >= 2) hs.add_node_list(s);
			sp1.clear();
		}

		qname = h.qname;
		hi = h.hi;

		if((h.flag & 0x4) >= 1) continue;

		vector<int64_t> v;
		h.get_matched_intervals(v);

		vector<int> sp2;
		for(int k = 0; k < v.size(); k++)
		{
			int32_t p1 = high32(v[k]);
			int32_t p2 = low32(v[k]);

			int k1 = locate_left_partial_exon(p1);
			int k2 = locate_right_partial_exon(p2);
			if(k1 < 0 || k2 < 0) continue;

			for(int j = k1; j <= k2; j++) sp2.push_back(j);
		}

		if(sp1.size() <= 0 || sp2.size() <= 0)
		{
			sp1.insert(sp1.end(), sp2.begin(), sp2.end());
			continue;
		}

		/*
		printf("sp2 = ( ");
		printv(sp2);
		printf(")\n");
		*/

		int x1 = -1, x2 = -1;
		if(h.isize < 0) 
		{
			x1 = sp1[max_element(sp1)];
			x2 = sp2[min_element(sp2)];
		}
		else
		{
			x1 = sp2[max_element(sp2)];
			x2 = sp1[min_element(sp1)];
		}

		vector<int> sp3;
		bool c = bridge_read(x1, x2, sp3);

		//printf("=========\n");

		if(c == false)
		{
			set<int> s(sp1.begin(), sp1.end());
			if(s.size() >= 2) hs.add_node_list(s);
			sp1 = sp2;
		}
		else
		{
			sp1.insert(sp1.end(), sp2.begin(), sp2.end());
			sp1.insert(sp1.end(), sp3.begin(), sp3.end());
		}
	}

	return 0;
}

bool bundle::bridge_read(int x, int y, vector<int> &v)
{
	v.clear();
	if(x >= y) return true;

	PEB e = gr.edge(x + 1, y + 1);
	if(e.second == true) return true;
	else return false;

	if(y - x >= 10) return false;

	long max = 9999999999;
	vector<long> table;
	vector<int> trace;
	int n = y - x + 1;
	table.resize(n, 0);
	trace.resize(n, -1);
	table[0] = 1;
	trace[0] = -1;
	for(int i = x + 1; i <= y; i++)
	{
		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.in_edges(i + 1); it1 != it2; it1++)
		{
			int s = (*it1)->source() - 1;
			int t = (*it1)->target() - 1;
			assert(t == i);
			if(s < x) continue;
			if(table[s - x] <= 0) continue;
			table[t - x] += table[s - x];
			trace[t - x] = s - x;
			if(table[t - x] >= max) return false;
		}
	}

	//printf("x = %d, y = %d, num-paths = %ld\n", x, y, table[n - 1]);
	if(table[n - 1] != 1) return false;

	//printf("path = ");

	v.clear();
	int p = n - 1;
	while(p >= 0)
	{
		p = trace[p];
		if(p <= 0) break;
		v.push_back(p + x);
		//printf("%d ", p + x);
	}
	//printf("\n");
	//assert(v.size() >= 1);

	return true;
}

int bundle::link_partial_exons()
{
	if(pexons.size() == 0) return 0;

	MPI lm;
	MPI rm;
	for(int i = 0; i < pexons.size(); i++)
	{
		int32_t l = pexons[i].lpos;
		int32_t r = pexons[i].rpos;

		assert(lm.find(l) == lm.end());
		assert(rm.find(r) == rm.end());
		lm.insert(PPI(l, i));
		rm.insert(PPI(r, i));
	}

	for(int i = 0; i < junctions.size(); i++)
	{
		junction &b = junctions[i];
		MPI::iterator li = rm.find(b.lpos);
		MPI::iterator ri = lm.find(b.rpos);

		assert(li != rm.end());
		assert(ri != lm.end());

		if(li != rm.end() && ri != lm.end())
		{
			b.lexon = li->second;
			b.rexon = ri->second;
		}
		else
		{
			b.lexon = b.rexon = -1;
		}
	}
	return 0;
}

int bundle::build_splice_graph()
{
	gr.clear();

	// vertices: start, each region, end
	gr.add_vertex();
	vertex_info vi0;
	vi0.lpos = lpos;
	vi0.rpos = lpos;
	gr.set_vertex_weight(0, 0);
	gr.set_vertex_info(0, vi0);
	for(int i = 0; i < pexons.size(); i++)
	{
		const partial_exon &r = pexons[i];
		int length = r.rpos - r.lpos;
		assert(length >= 1);
		gr.add_vertex();
		gr.set_vertex_weight(i + 1, r.ave < 1.0 ? 1.0 : r.ave);
		vertex_info vi;
		vi.lpos = r.lpos;
		vi.rpos = r.rpos;
		vi.length = length;
		vi.stddev = r.dev;// < 1.0 ? 1.0 : r.dev;
		//vi.adjust = r.adjust;
		gr.set_vertex_info(i + 1, vi);
	}

	gr.add_vertex();
	vertex_info vin;
	vin.lpos = rpos;
	vin.rpos = rpos;
	gr.set_vertex_weight(pexons.size() + 1, 0);
	gr.set_vertex_info(pexons.size() + 1, vin);

	// edges: each junction => and e2w
	for(int i = 0; i < junctions.size(); i++)
	{
		const junction &b = junctions[i];

		if(b.lexon < 0 || b.rexon < 0) continue;

		const partial_exon &x = pexons[b.lexon];
		const partial_exon &y = pexons[b.rexon];

		edge_descriptor p = gr.add_edge(b.lexon + 1, b.rexon + 1);
		assert(b.count >= 1);
		edge_info ei;
		ei.weight = b.count;
		gr.set_edge_info(p, ei);
		gr.set_edge_weight(p, b.count);
	}

	// edges: connecting start/end and pexons
	int ss = 0;
	int tt = pexons.size() + 1;
	for(int i = 0; i < pexons.size(); i++)
	{
		const partial_exon &r = pexons[i];

		if(r.ltype == START_BOUNDARY)
		{
			edge_descriptor p = gr.add_edge(ss, i + 1);
			double w = r.ave;
			if(i >= 1 && pexons[i - 1].rpos == r.lpos) w -= pexons[i - 1].ave;
			if(w < 1.0) w = 1.0;
			gr.set_edge_weight(p, w);
			edge_info ei;
			ei.weight = w;
			gr.set_edge_info(p, ei);
		}

		if(r.rtype == END_BOUNDARY) 
		{
			edge_descriptor p = gr.add_edge(i + 1, tt);
			double w = r.ave;
			if(i < pexons.size() - 1 && pexons[i + 1].lpos == r.rpos) w -= pexons[i + 1].ave;
			if(w < 1.0) w = 1.0;
			gr.set_edge_weight(p, w);
			edge_info ei;
			ei.weight = w;
			gr.set_edge_info(p, ei);
		}
	}

	// edges: connecting adjacent pexons => e2w
	for(int i = 0; i < (int)(pexons.size()) - 1; i++)
	{
		const partial_exon &x = pexons[i];
		const partial_exon &y = pexons[i + 1];

		if(x.rpos != y.lpos) continue;

		assert(x.rpos == y.lpos);
		
		int xd = gr.out_degree(i + 1);
		int yd = gr.in_degree(i + 2);
		double wt = (xd < yd) ? x.ave : y.ave;
		//int32_t xr = compute_overlap(mmap, x.rpos - 1);
		//int32_t yl = compute_overlap(mmap, y.lpos);
		//double wt = xr < yl ? xr : yl;

		edge_descriptor p = gr.add_edge(i + 1, i + 2);
		double w = (wt < 1.0) ? 1.0 : wt;
		gr.set_edge_weight(p, w);
		edge_info ei;
		ei.weight = w;
		gr.set_edge_info(p, ei);
	}


	return 0;
}

int bundle::extend_isolated_end_boundaries()
{
	for(int i = 1; i < gr.num_vertices(); i++)
	{
		if(gr.in_degree(i) != 1) continue;
		if(gr.out_degree(i) != 1) continue;

		edge_iterator it1, it2;
		tie(it1, it2) = gr.in_edges(i);
		edge_descriptor e1 = (*it1);
		tie(it1, it2) = gr.out_edges(i);
		edge_descriptor e2 = (*it1);
		int s = e1->source();
		int t = e2->target();

		if(gr.out_degree(s) != 1) continue;
		if(t != gr.num_vertices() - 1) continue;
		if(gr.get_edge_weight(e1) >= 2.5) continue;
		if(gr.get_edge_weight(e2) >= 5.0) continue;
		double w1 = gr.get_edge_weight(e1);
		double wv = gr.get_vertex_weight(s);
		if(wv < w1 * w1 + 5) continue;

		//if(gr.get_vertex_weight(s) <= 5.0) continue;
		if(gr.get_vertex_info(s).rpos == gr.get_vertex_info(i).lpos) continue;

		double w = gr.get_vertex_weight(s) - gr.get_edge_weight(e1);
		edge_descriptor e = gr.add_edge(s, t);
		gr.set_edge_weight(e, w);
		gr.set_edge_info(e, edge_info());
		
		printf("extend isolated end boundary: (%d, %.2lf) -- (%.2lf) -- (%d, %.2lf)\n", s, gr.get_vertex_weight(s),
				gr.get_edge_weight(e1), i, gr.get_vertex_weight(i));

	}
	return 0;
}

int bundle::extend_isolated_start_boundaries()
{
	for(int i = 1; i < gr.num_vertices(); i++)
	{
		if(gr.in_degree(i) != 1) continue;
		if(gr.out_degree(i) != 1) continue;

		edge_iterator it1, it2;
		tie(it1, it2) = gr.in_edges(i);
		edge_descriptor e1 = (*it1);
		tie(it1, it2) = gr.out_edges(i);
		edge_descriptor e2 = (*it1);
		int s = e1->source();
		int t = e2->target();

		if(s != 0) continue;
		if(gr.in_degree(t) != 1) continue;
		if(gr.get_edge_weight(e2) >= 2.5) continue;
		if(gr.get_edge_weight(e1) >= 5.0) continue;
		double w2 = gr.get_edge_weight(e2);
		double wv = gr.get_vertex_weight(t);
		if(wv < w2 * w2 + 5) continue;

		if(gr.get_vertex_info(i).rpos == gr.get_vertex_info(t).lpos) continue;

		double w = gr.get_vertex_weight(t) - gr.get_edge_weight(e2);
		edge_descriptor e = gr.add_edge(s, t);
		gr.set_edge_weight(e, w);
		gr.set_edge_info(e, edge_info());
		
		printf("extend isolated start boundary: (%d, %.2lf) -- (%.2lf) -- (%d, %.2lf)\n", i, gr.get_vertex_weight(i),
				gr.get_edge_weight(e2), t, gr.get_vertex_weight(t));

		gr.remove_edge(e1);
		gr.remove_edge(e2);
	}
	return 0;
}

int bundle::remove_inner_start_boundaries()
{
	for(int i = 1; i < gr.num_vertices(); i++)
	{
		if(gr.in_degree(i) != 1) continue;
		if(gr.out_degree(i) != 1) continue;

		edge_iterator it1, it2;
		tie(it1, it2) = gr.in_edges(i);
		edge_descriptor e1 = (*it1);
		tie(it1, it2) = gr.out_edges(i);
		edge_descriptor e2 = (*it1);
		int s = e1->source();
		int t = e2->target();
		vertex_info vi = gr.get_vertex_info(i);
		double wv = gr.get_vertex_weight(i);
		double ww = gr.get_vertex_weight(t);

		if(s != 0) continue;
		if(gr.in_degree(t) == 1) continue;

		bool b1 = true;
		if(vi.stddev >= 0.01) b1 = false;

		bool b2 = true;
		if(wv > min_inner_boundary_weight) b2 = false;
		if(t != i + 1) b2 = false;
		if(vi.rpos != gr.get_vertex_info(t).lpos) b2 = false;
		//if(ww < 1.5 * wv) continue;

		if(b1 == false && b2 == false) continue;

		printf("remove inner start boundary: vertex = %d, weight = %.2lf, length = %d\n", i, wv, vi.length);

		gr.clear_vertex(i);
	}
	return 0;
}

int bundle::remove_inner_end_boundaries()
{
	for(int i = 1; i < gr.num_vertices(); i++)
	{
		if(gr.in_degree(i) != 1) continue;
		if(gr.out_degree(i) != 1) continue;

		edge_iterator it1, it2;
		tie(it1, it2) = gr.in_edges(i);
		edge_descriptor e1 = (*it1);
		tie(it1, it2) = gr.out_edges(i);
		edge_descriptor e2 = (*it1);
		int s = e1->source();
		int t = e2->target();
		vertex_info vi = gr.get_vertex_info(i);
		double wv = gr.get_vertex_weight(i);
		double ww = gr.get_vertex_weight(s);

		if(t != gr.num_vertices() - 1) continue;
		if(gr.out_degree(s) == 1) continue;

		//if(ww < 1.5 * wv) continue;
	
		bool b1 = true;
		if(vi.stddev >= 0.01) b1 = false;

		bool b2 = true;
		if(vi.lpos != gr.get_vertex_info(s).rpos) b2 = false;
		if(i != s + 1) b2 = false;
		if(wv > min_inner_boundary_weight) b2 = false;
	
		if(b1 == false && b2 == false) continue;

		printf("remove inner end boundary: vertex = %d, weight = %.2lf, length = %d\n", i, wv, vi.length);

		gr.clear_vertex(i);
	}
	return 0;
}

int bundle::remove_inner_vertices()
{
	for(int i = 1; i < gr.num_vertices(); i++)
	{
		if(gr.in_degree(i) != 1) continue;
		if(gr.out_degree(i) != 1) continue;

		edge_iterator it1, it2;
		tie(it1, it2) = gr.in_edges(i);
		edge_descriptor e1 = (*it1);
		tie(it1, it2) = gr.out_edges(i);
		edge_descriptor e2 = (*it1);
		int s = e1->source();
		int t = e2->target();
		double wv = gr.get_vertex_weight(i);
		vertex_info vi = gr.get_vertex_info(i);

		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;
		if(gr.get_vertex_info(s).rpos != vi.lpos) continue;
		if(gr.get_vertex_info(t).lpos != vi.rpos) continue;

		PEB p = gr.edge(s, t);
		if(p.second == false) continue;

		edge_descriptor ee = p.first;
		double we = gr.get_edge_weight(ee);

		if(wv > we) continue;
		if(wv > min_inner_vertex_weight) continue;

		printf("clear inner exon %d, weight = %.2lf, length = %d, edge weight = %.2lf\n", i, wv, vi.length, we);

		gr.clear_vertex(i);
	}
	return 0;
}

int bundle::join_single_exon_transcripts()
{
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		int x = i;
		int y = i + 1;
		if(gr.in_degree(x) != 1) continue;
		if(gr.in_degree(y) != 1) continue;
		if(gr.out_degree(x) != 1) continue;
		if(gr.out_degree(y) != 1) continue;

		edge_iterator it1, it2;

		tie(it1, it2) = gr.in_edges(x);
		edge_descriptor e1 = (*it1);
		tie(it1, it2) = gr.out_edges(x);
		edge_descriptor e2 = (*it1);

		tie(it1, it2) = gr.in_edges(y);
		edge_descriptor e3 = (*it1);
		tie(it1, it2) = gr.out_edges(y);
		edge_descriptor e4 = (*it1);

		if(e1->source() != 0) continue;
		if(e2->target() != gr.num_vertices() - 1) continue;
		if(e3->source() != 0) continue;
		if(e4->target() != gr.num_vertices() - 1) continue;

		int32_t p1 = gr.get_vertex_info(x).rpos;
		int32_t p2 = gr.get_vertex_info(y).lpos;

		if(p1 >= p2 || p2 - p1 > min_bundle_gap) continue;

		double w2 = gr.get_edge_weight(e2);
		double w3 = gr.get_edge_weight(e3);

		gr.remove_edge(e2);
		gr.remove_edge(e3);

		edge_descriptor ee = gr.add_edge(x, y);
		double ww = 0.5 * w2 + 0.5 * w3;
		gr.set_edge_weight(ee, ww);
		gr.set_edge_info(ee, edge_info());
	}
	return 0;
}

VE bundle::compute_maximal_edges()
{
	typedef pair<double, edge_descriptor> PDE;
	vector<PDE> ve;

	undirected_graph ug;
	edge_iterator it1, it2;
	for(int i = 0; i < gr.num_vertices(); i++) ug.add_vertex();
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		double w = gr.get_edge_weight(e);
		int s = e->source();
		int t = e->target();
		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;
		ug.add_edge(s, t);
		ve.push_back(PDE(w, e));
	}

	vector<int> vv = ug.assign_connected_components();

	sort(ve.begin(), ve.end());

	for(int i = 1; i < ve.size(); i++) assert(ve[i - 1].first <= ve[i].first);

	VE x;
	set<int> sc;
	for(int i = ve.size() - 1; i >= 0; i--)
	{
		edge_descriptor e = ve[i].second;
		double w = gr.get_edge_weight(e);
		if(w < 1.5) break;
		int s = e->source();
		int t = e->target();
		if(s == 0) continue;
		if(t == gr.num_vertices()) continue;
		int c1 = vv[s];
		int c2 = vv[t];
		assert(c1 == c2);
		if(sc.find(c1) != sc.end()) continue;
		x.push_back(e);
		sc.insert(c1);
	}
	return x;
}

int bundle::refine_splice_graph()
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

int bundle::remove_small_edges()
{
	set<int> sv1;
	set<int> sv2;
	SE se;
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		double w = gr.get_edge_weight(*it1);
		int32_t p1 = gr.get_vertex_info(s).rpos;
		int32_t p2 = gr.get_vertex_info(t).lpos;
		if(w < min_splice_edge_weight) continue;
		//if(p1 == p2 && w < min_consecutive_edge_weight) continue;
		se.insert(*it1);
		sv1.insert(t);
		sv2.insert(s);
	}

	VE me = compute_maximal_edges();
	for(int i = 0; i < me.size(); i++)
	{
		edge_descriptor ee = me[i];
		se.insert(ee);
		sv1.insert(ee->target());
		sv2.insert(ee->source());
	}

	while(true)
	{
		bool b = false;
		for(SE::iterator it = se.begin(); it != se.end(); it++)
		{
			edge_descriptor e = (*it);
			int s = e->source(); 
			int t = e->target();
			if(sv1.find(s) == sv1.end() && s != 0)
			{
				edge_descriptor ee = gr.max_in_edge(s);
				assert(ee != null_edge);
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(s);
				sv2.insert(ee->source());
				b = true;
			}
			if(sv2.find(t) == sv2.end() && t != gr.num_vertices() - 1)
			{
				edge_descriptor ee = gr.max_out_edge(t);
				assert(ee != null_edge);
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(ee->target());
				sv2.insert(t);
				b = true;
			}
			if(b == true) break;
		}
		if(b == false) break;
	}

	VE ve;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		if(se.find(*it1) != se.end()) continue;
		ve.push_back(*it1);
	}

	for(int i = 0; i < ve.size(); i++)
	{
		printf("remove edge (%d, %d), weight = %.2lf\n", ve[i]->source(), ve[i]->target(), gr.get_edge_weight(ve[i]));
		gr.remove_edge(ve[i]);
	}

	return 0;
}

int bundle::count_junctions() const
{
	int x = 0;
	for(int i = 0; i < junctions.size(); i++)
	{
		x += junctions[i].count;
	}
	return x;
}

int bundle::print(int index)
{
	printf("\nBundle %d: ", index);

	// statistic xs
	int n0 = 0, np = 0, nq = 0;
	for(int i = 0; i < hits.size(); i++)
	{
		if(hits[i].xs == '.') n0++;
		if(hits[i].xs == '+') np++;
		if(hits[i].xs == '-') nq++;
	}

	printf("tid = %d, #hits = %lu, #partial-exons = %lu, range = %s:%d-%d, orient = %c (%d, %d, %d)\n",
			tid, hits.size(), pexons.size(), chrm.c_str(), lpos, rpos, strand, n0, np, nq);

	// print hits
	//for(int i = 0; i < hits.size(); i++) hits[i].print();

	// print regions
	for(int i = 0; i < regions.size(); i++)
	{
		regions[i].print(i);
	}

	// print junctions 
	for(int i = 0; i < junctions.size(); i++)
	{
		junctions[i].print(i);
	}

	// print partial exons
	for(int i = 0; i < pexons.size(); i++)
	{
		pexons[i].print(i);
	}

	// print hyper-edges
	//hs.print();

	return 0;
}

int bundle::output_transcripts(ofstream &fout, const vector<path> &p, const string &gid) const
{
	for(int i = 0; i < p.size(); i++)
	{
		string tid = gid + "." + tostring(i);
		output_transcript(fout, p[i], gid, tid);
	}
	return 0;
}

int bundle::output_transcript(ofstream &fout, const path &p, const string &gid, const string &tid) const
{
	fout.precision(2);
	fout<<fixed;

	const vector<int> &v = p.v;
	double coverage = p.abd;		// number of molecular
	double numreads = p.reads / average_read_length;

	assert(v[0] == 0);
	assert(v[v.size() - 1] == pexons.size() + 1);
	if(v.size() < 2) return 0;

	int ss = v[1];
	int tt = v[v.size() - 2];
	int32_t ll = pexons[ss - 1].lpos;
	int32_t rr = pexons[tt - 1].rpos;

	fout<<chrm.c_str()<<"\t";		// chromosome name
	fout<<algo.c_str()<<"\t";		// source
	fout<<"transcript\t";			// feature
	fout<<ll + 1<<"\t";				// left position
	fout<<rr<<"\t";					// right position
	fout<<1000<<"\t";				// score, now as abundance
	fout<<strand<<"\t";				// strand
	fout<<".\t";					// frame
	fout<<"gene_id \""<<gid.c_str()<<"\"; ";
	fout<<"transcript_id \""<<tid.c_str()<<"\"; ";
	fout<<"coverage \""<<coverage<<"\"; ";
	fout<<"numreads \""<<numreads<<"\";"<<endl;

	join_interval_map jmap;
	for(int k = 1; k < v.size() - 1; k++)
	{
		const partial_exon &r = pexons[v[k] - 1];
		jmap += make_pair(ROI(r.lpos, r.rpos), 1);
	}

	int cnt = 0;
	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		fout<<chrm.c_str()<<"\t";			// chromosome name
		fout<<algo.c_str()<<"\t";			// source
		fout<<"exon\t";						// feature
		fout<<lower(it->first) + 1<<"\t";	// left position
		fout<<upper(it->first)<<"\t";		// right position
		fout<<1000<<"\t";					// score
		fout<<strand<<"\t";					// strand
		fout<<".\t";						// frame
		fout<<"gene_id \""<<gid.c_str()<<"\"; ";
		fout<<"transcript_id \""<<tid.c_str()<<"\"; ";
		fout<<"exon_number \""<<++cnt<<"\"; ";
		fout<<"coverage \""<<coverage<<"\"; ";
		fout<<"numreads \""<<numreads<<"\";"<<endl;
	}
	return 0;
}

int bundle::output_transcripts(gene &gn, const vector<path> &p, const string &gid) const
{
	for(int i = 0; i < p.size(); i++)
	{
		string tid = gid + "." + tostring(i);
		transcript trst;
		output_transcript(trst, p[i], gid, tid);
		gn.add_transcript(trst);
	}
	return 0;
}

int bundle::output_transcript(transcript &trst, const path &p, const string &gid, const string &tid) const
{
	trst.seqname = chrm;
	trst.source = algo;
	trst.gene_id = gid;
	trst.transcript_id = tid;
	trst.coverage = p.abd;
	trst.numreads = p.reads / average_read_length;
	trst.strand = strand;

	const vector<int> &v = p.v;
	join_interval_map jmap;
	for(int k = 1; k < v.size() - 1; k++)
	{
		const partial_exon &r = pexons[v[k] - 1];
		jmap += make_pair(ROI(r.lpos, r.rpos), 1);
	}
	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		trst.add_exon(lower(it->first), upper(it->first));
	}
	return 0;
}
