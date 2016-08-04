#include <cassert>
#include <cstdio>
#include <map>
#include <iomanip>
#include <fstream>

#include "bundle.h"
#include "binomial.h"
#include "region.h"

bundle::bundle(const bundle_base &bb)
	: bundle_base(bb)
{
}

bundle::~bundle()
{}

int bundle::build()
{
	check_left_ascending();
	infer_junctions();
	build_junction_graph();
	//infer_hyper_junctions();
	//draw_junction_graph("jr.tex");
	//test_junction_graph();
	process_hits();
	build_super_regions();
	build_partial_exons();
	link_partial_exons();
	infer_hyper_edges();
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

int bundle::infer_junctions()
{
	map<int64_t, junction> m;
	vector<int64_t> v;
	for(int i = 0; i < hits.size(); i++)
	{
		hits[i].get_splice_positions(v);
		if(v.size() == 0) continue;
		for(int k = 0; k < v.size(); k++)
		{
			int64_t p = v[k];
			if(m.find(p) == m.end()) 
			{
				junction sp(p, 1, hits[i].qual, hits[i].qual);
				m.insert(pair<int64_t, junction>(p, sp));
			}
			else
			{
				m[p].count++;
				if(hits[i].qual < m[p].min_qual) m[p].min_qual = hits[i].qual;
				if(hits[i].qual > m[p].max_qual) m[p].max_qual = hits[i].qual;
			}
		}
	}

	map<int64_t, junction>::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		if(it->second.count < min_splice_boundary_hits) continue;
		if(it->second.max_qual < min_max_splice_boundary_qual) continue;
		junctions.push_back(it->second);
	}
	return 0;
}

int bundle::infer_hyper_junctions()
{
	map<int64_t, int> p2i;
	for(int i = 0; i < junctions.size(); i++)
	{
		junction &jc = junctions[i];
		int64_t p = pack(jc.lpos, jc.rpos);
		p2i.insert(pair<int64_t, int>(p, i));
	}

	vector<int> list;
	map<string, PI> m;
	vector<int64_t> v;
	for(int i = 0; i < hits.size(); i++)
	{
		hits[i].get_splice_positions(v);
		if(v.size() == 0) continue;
		string s = hits[i].qname;

		if(m.find(s) == m.end())
		{
			PI pi(list.size(), list.size());
			for(int k = 0; k < v.size(); k++)
			{
				int64_t p = v[k];
				if(p2i.find(p) == p2i.end()) continue;
				list.push_back(p2i[p]);
				pi.second++;
			}
			if(pi.second > pi.first)
			{
				m.insert(pair<string, PI>(s, pi));
			}
		}
		else
		{
			int l = m[s].first;
			int r = m[s].second;
			vector<int> vv1(list.begin() + l, list.begin() + r);
			vector<int> vv2;
			for(int k = 0; k < v.size(); k++)
			{
				int64_t p = v[k];
				if(p2i.find(p) == p2i.end()) continue;
				vv2.push_back(p2i[p]);

				bool b = true;
				for(int j = l; j < r; j++)
				{
					if(list[j] == p2i[p]) b = false;
					if(b == false) break;
				}
				//if(b == true) vv1.push_back(p2i[p]);
			}

			// super junctions
			printf("super junction: ");
			printv(vv1);
			printf(" | ");
			printv(vv2);
			printf("\n");
		}
	}
	return 0;
}

int bundle::build_junction_graph()
{
	map<int, int> s;
	s.insert(PI(lpos, START_BOUNDARY));
	s.insert(PI(rpos, END_BOUNDARY));
	for(int i = 0; i < junctions.size(); i++)
	{
		int l = junctions[i].lpos;
		int r = junctions[i].rpos;
		if(s.find(l) == s.end()) s.insert(PI(l, LEFT_SPLICE));
		if(s.find(r) == s.end()) s.insert(PI(r, RIGHT_SPLICE));
	}
	vector<PI> v(s.begin(), s.end());
	sort(v.begin(), v.end());

	// position to vertex
	map<int, int> p2v;
	for(int i = 0; i < v.size(); i++)
	{
		p2v.insert(PI(v[i].first, i));
	}

	// vertices
	jr.clear();
	for(int i = 0; i < v.size(); i++)
	{
		jr.add_vertex();
		vertex_info vi;
		vi.type = v[i].second;
		jr.set_vertex_info(i, vi);
		jr.set_vertex_weight(i, v[i].first);
	}

	// edges: connecting adjacent regions
	for(int i = 0; i < v.size() - 1; i++)
	{
		edge_descriptor p = jr.add_edge(i, i + 1);
		jr.set_edge_weight(p, -1);
	}

	// edges: each junction
	for(int i = 0; i < junctions.size(); i++)
	{
		int l = junctions[i].lpos;
		int r = junctions[i].rpos;
		assert(p2v.find(l) != p2v.end());
		assert(p2v.find(r) != p2v.end());
		edge_descriptor p = jr.add_edge(p2v[l], p2v[r]);
		jr.set_edge_weight(p, i);
	}

	return 0;
}

int bundle::draw_junction_graph(const string &file)
{
	char buf[10240];

	MIS mis;
	for(int i = 0; i < jr.num_vertices(); i++)
	{
		double w = jr.get_vertex_weight(i);
		sprintf(buf, "%.0lf", w);
		mis.insert(PIS(i, buf));
	}

	MES mes;
	edge_iterator it1, it2;
	for(tie(it1, it2) = jr.edges(); it1 != it2; it1++)
	{
		double w = jr.get_edge_weight(*it1);
		sprintf(buf, "%.0lf", w);
		mes.insert(PES(*it1, buf));
	}

	jr.draw(file, mis, mes, 3.0);
	return 0;
}

int bundle::test_junction_graph()
{
	for(int i = 0; i < jr.num_vertices(); i++)
	{
		for(int j = i + 1; j < jr.num_vertices(); j++)
		{
			VE ve;
			int p = traverse_junction_graph1(i, j, ve);
			printf("path1 from %d to %d, length = %d :", i, j, p);
			for(int k = 0; k < ve.size(); k++)
			{
				int s = ve[k]->source();
				int t = ve[k]->target();
				int id = (int)(jr.get_edge_weight(ve[k]));
				printf("(%d, %d, %d) <- ", s, t, id);
			}
			printf(" START\n");
		}
	}
	return 0;
}

int bundle::search_junction_graph(int32_t p)
{
	assert(jr.num_vertices() >= 2);
	int l = 0;
	int r = jr.num_vertices() - 1;
	while(l < r)
	{
		int m = (l + r) / 2;
		int32_t p1 = (int32_t)(jr.get_vertex_weight(m));
		int32_t p2 = (int32_t)(jr.get_vertex_weight(m + 1));
		assert(p1 < p2);
		if(p >= p1 && p < p2) return m;
		if(p < p1) r = m;
		if(p >= p2) l = m + 1;
	}
	return -1;
}

int bundle::traverse_junction_graph1(int s, int t)
{
	VE ve;
	return traverse_junction_graph1(s, t, ve);
}

int bundle::traverse_junction_graph1(int s, int t, VE &ve)
{
	ve.clear();
	if(s > t) return -1;

	vector<double> v1, v2;		// shortest path/with one edge
	VE ve1, ve2;				// tracing back pointers
	v1.push_back(0);
	v2.push_back(-1);
	ve1.push_back(null_edge);
	ve2.push_back(null_edge);

	for(int k = s + 1; k <= t; k++)
	{
		edge_descriptor ee1 = null_edge;
		edge_descriptor ee2 = null_edge;
		double ww1 = DBL_MAX;
		double ww2 = DBL_MAX;
		edge_iterator it1, it2;
		for(tie(it1, it2) = jr.in_edges(k); it1 != it2; it1++)
		{
			int ss = (*it1)->source();
			if(ss < s) continue;

			double w = 0;
			int id = (int)(jr.get_edge_weight(*it1));
			if(id < 0) 
			{
				assert(ss == k - 1);
				w = jr.get_vertex_weight(k) - jr.get_vertex_weight(ss);

				if(v1[ss - s] + w < ww1)
				{
					ww1 = v1[ss - s] + w;
					ee1 = *it1;
				}

				if(v1[ss - s] + w < ww2)
				{
					ww2 = v1[ss - s] + w;
					ee2 = *it1;
				}
			}
			else
			{
				if(v1[ss - s] < ww1)
				{
					ww1 = v1[ss - s];
					ee1 = *it1;
				}

				if(v2[ss - s] >= 0 && v2[ss - s] < ww2)
				{
					ww2 = v2[ss - s];
					ee2 = *it1;
				}
			}
		}

		assert(ee1 != null_edge);
		v1.push_back(ww1);
		ve1.push_back(ee1);

		if(ee2 == null_edge)
		{
			v2.push_back(-1);
			ve2.push_back(null_edge);
		}
		else
		{
			v2.push_back(ww2);
			ve2.push_back(ee2);
		}
	}

	if(v2[t - s] <= 0) return -1;

	int k = t - s;
	while(ve2[k] != null_edge)
	{
		ve.push_back(ve2[k]);
		int id = (int)(jr.get_edge_weight(ve2[k]));
		k = ve2[k]->source() - s;
		if(id < 0) break;
	}

	while(ve1[k] != null_edge)
	{
		ve.push_back(ve1[k]);
		k = ve1[k]->source() - s;
	}

	return (int)(v2[t - s]);
}

int bundle::traverse_junction_graph(int s, int t, VE &ve)
{
	ve.clear();
	if(s > t) return -1;
	vector<double> v1;
	VE v2;
	v1.push_back(0);
	v2.push_back(null_edge);

	for(int k = s + 1; k <= t; k++)
	{
		edge_descriptor ee = null_edge;
		double ww = DBL_MAX;
		edge_iterator it1, it2;
		for(tie(it1, it2) = jr.in_edges(k); it1 != it2; it1++)
		{
			int ss = (*it1)->source();
			if(ss < s) continue;

			double w = 0;
			int id = (int)(jr.get_edge_weight(*it1));
			if(id < 0) 
			{
				assert(ss == k - 1);
				w = jr.get_vertex_weight(k) - jr.get_vertex_weight(ss);
			}

			if(v1[ss - s] + w < ww)
			{
				ww = v1[ss - s] + w;
				ee =  *it1;
			}
		}

		assert(ee != null_edge);
		v1.push_back(ww);
		v2.push_back(ee);
	}

	int k = t - s;
	while(v2[k] != null_edge)
	{
		ve.push_back(v2[k]);
		k = v2[k]->source() - s;
	}
	return (int)(v1[t - s]);
}

int bundle::process_hits()
{
	imap.clear();
	for(int i = 0; i < hits.size(); i++)
	{
		hit &h = hits[i];

		if(use_paired_end == false)
		{
			add_mapped_intervals(h, h.rpos);
			continue;
		}

		if(h.isize <= 0)
		{
			// the second segment
			add_mapped_intervals(h, h.rpos);
		}
		else if(h.rpos > h.mpos)
		{
			// two segments overlap
			add_mapped_intervals(h, h.mpos);
		}
		else
		{
			// two independent segments
			add_mapped_intervals(h, h.rpos);
			add_gapped_intervals(h);
		}
	}
	return 0;
}

int bundle::add_mapped_intervals(const hit &h, int32_t rr)
{
	vector<int64_t> v;
	h.get_matched_intervals(v);
	for(int k = 0; k < v.size(); k++)
	{
		int32_t s = high32(v[k]);
		int32_t t = low32(v[k]);
		if(s >= rr) s = rr;
		if(t >= rr) t = rr;
		imap += make_pair(ROI(s, t), 1);
	}
	return 0;
}

int bundle::add_gapped_intervals(const hit &h)
{
	assert(h.rpos <= h.mpos);
	int li = search_junction_graph(h.rpos - 1);
	if(li < 0) return -1;
	assert(li >= 0 && li < jr.num_vertices() - 1);
	int32_t l1 = (int32_t)(jr.get_vertex_weight(li));
	int32_t l2 = (int32_t)(jr.get_vertex_weight(li + 1));
	assert(h.rpos - 1 >= l1 && h.rpos - 1 < l2);
	
	if(h.mpos < l2)
	{
		imap += make_pair(ROI(h.rpos, h.mpos), 1);
		return 0;
	}

	int ri = search_junction_graph(h.mpos);
	if(ri < 0) return -1;
	assert(ri >= 0 && ri < jr.num_vertices() - 1);
	assert(ri > li);
	int32_t r1 = (int32_t)(jr.get_vertex_weight(ri));
	int32_t r2 = (int32_t)(jr.get_vertex_weight(ri + 1));
	assert(h.mpos >= r1 && h.mpos < r2);

	imap += make_pair(ROI(h.rpos, l2), 1);
	imap += make_pair(ROI(r1, h.mpos), 1);

	VE ve;
	int ss = traverse_junction_graph(li + 1, ri, ve);

	//int ms = (int)(ave_isize) - (l2 - h.rpos) - (h.mpos - r1);
	//if(fabs(ms - ss) >= 100) return 0;
	//printf("span hit: remaining insert size = %d, best insert size = %d, li = %d, ri = %d\n", ms, ss, li, ri);

	assert(ss >= 0);

	for(int k = 0; k < ve.size(); k++)
	{
		edge_descriptor e = ve[k];
		int id = jr.get_edge_weight(e);
		if(id < 0)
		{
			int32_t l = (int32_t)jr.get_vertex_weight(e->source());
			int32_t r = (int32_t)jr.get_vertex_weight(e->target());
			imap += make_pair(ROI(l, r), 1);
		}
		else
		{
			assert(id >= 0 && id < junctions.size());
			junctions[id].count++;
		}
	}
	
	return 0;
}

int bundle::build_super_region(int s, super_region &sr)
{
	assert(s >= 0 && s < jr.num_vertices() - 1);

	// find the right position
	int t = s + 1;
	edge_iterator it1, it2;
	for(; t < jr.num_vertices() - 1; t++)
	{
		bool b = true;
		for(tie(it1, it2) = jr.in_edges(t); it1 != it2; it1++)
		{
			if((*it1)->source() != t - 1) b = false;
			if(b == false) break;
		}
		if(b == false) break;

		for(tie(it1, it2) = jr.out_edges(t); it1 != it2; it1++)
		{
			if((*it1)->target() != t + 1) b = false;
			if(b == false) break;
		}
		if(b == false) break;
	}

	sr.clear();

	if((t - s) % 2 == 0) t--;

	for(int i = s; i < t; i += 2)
	{
		// check region (i, i + 1)
		int32_t lpos = (int32_t)(jr.get_vertex_weight(i));
		int32_t rpos = (int32_t)(jr.get_vertex_weight(i + 1));
		int ltype = jr.get_vertex_info(i).type; 
		int rtype = jr.get_vertex_info(i + 1).type; 

		region r(lpos, rpos, ltype, rtype, &imap);

		if(r.empty == true) return 0;

		sr.add_region(r);

		if(i + 1 >= t) return 0;
		if(i != s && jr.out_degree(i) != 1) return 0;
		if(i != t - 1 && jr.in_degree(i + 1) != 1) return 0;

		// check region (i + 1, i + 2)
		if(jr.out_degree(i + 1) != 2) return 0;
		if(jr.in_degree(i + 2) != 2) return 0;

		lpos = (int32_t)(jr.get_vertex_weight(i + 1));
		rpos = (int32_t)(jr.get_vertex_weight(i + 2));
		ltype = jr.get_vertex_info(i + 1).type; 
		rtype = jr.get_vertex_info(i + 2).type; 

		region r0(lpos, rpos, ltype, rtype, &imap);

		if(r0.empty == false) return 0;
	}
	return 0;
}

int bundle::build_super_regions()
{
	for(int k = 0; k < jr.num_vertices() - 1;)
	{
		super_region sr(&imap);
		build_super_region(k, sr);
		if(sr.size() >= 1) srs.push_back(sr);

		if(sr.size() == 0) k++;
		else k += sr.size() * 2 - 1;
	}
	return 0;
}

int bundle::build_partial_exons()
{
	for(int i = 0; i < srs.size(); i++)
	{
		super_region &sr = srs[i];
		sr.build();
		pexons.insert(pexons.end(), sr.pexons.begin(), sr.pexons.end());
	}
	return 0;
}

int bundle::search_partial_exons(int32_t x)
{
	int l = 0;
	int r = pexons.size() - 1;
	while(l <= r)
	{
		int m = (l + r) / 2;
		partial_exon &p = pexons[m];
		if(x >= p.lpos && x < p.rpos) return m;

		if(x < p.lpos) r = m - 1;
		if(x >= p.rpos) l = m + 1;
	}
	return -1;
}

int bundle::infer_hyper_edges()
{
	vector<int> list;
	map<string, PI> m;
	vector<int64_t> v;
	map<hyper_edge, int> mhe;
	for(int i = 0; i < hits.size(); i++)
	{
		hit &h = hits[i];
		h.get_matched_intervals(v);
		if(v.size() == 0) continue;
		string s = h.qname;

		if(h.isize > 0 && h.rpos < h.mpos)
		{
			int li = search_junction_graph(h.rpos - 1);
			if(li < 0) continue;
			assert(li >= 0 && li < jr.num_vertices() - 1);

			int32_t l1 = (int32_t)(jr.get_vertex_weight(li));
			int32_t l2 = (int32_t)(jr.get_vertex_weight(li + 1));
			assert(h.rpos - 1 >= l1 && h.rpos - 1 < l2);

			if(h.mpos > l2)
			{
				int ri = search_junction_graph(h.mpos);
				if(ri < 0) continue;

				assert(ri >= 0 && ri < jr.num_vertices() - 1);
				assert(ri > li);
				int32_t r1 = (int32_t)(jr.get_vertex_weight(ri));
				int32_t r2 = (int32_t)(jr.get_vertex_weight(ri + 1));
				assert(h.mpos >= r1 && h.mpos < r2);
				
				int d1 = r1 - l2;
				assert(d1 >= 0);

				int d2 = traverse_junction_graph1(li, ri);
				
				if(d2 >= 0 && d1 > d2) continue;
			}
		}

		if(m.find(s) == m.end())
		{
			PI pi(list.size(), list.size());
			for(int k = 0; k < v.size(); k++)
			{
				int32_t p1 = high32(v[k]);
				int32_t p2 = low32(v[k]);

				int e1 = search_partial_exons(p1);
				int e2 = search_partial_exons(p2 - 1);

				if(e1 <= -1 || e2 <= -1) continue;

				for(int j = e1; j <= e2; j++) 
				{
					list.push_back(j);
					pi.second++;
				}
				if(pi.second > pi.first)
				{
					m.insert(pair<string, PI>(s, pi));
				}
			}
		}
		else
		{
			int l = m[s].first;
			int r = m[s].second;
			set<int> sp(list.begin() + l, list.begin() + r);
			for(int k = 0; k < v.size(); k++)
			{
				int32_t p1 = high32(v[k]);
				int32_t p2 = low32(v[k]);

				int e1 = search_partial_exons(p1);
				int e2 = search_partial_exons(p2 - 1);

				if(e1 <= -1 || e2 <= -1) continue;

				for(int j = e1; j <= e2; j++) 
				{
					if(sp.find(j) == sp.end()) sp.insert(j);
				}
			}

			// super junctions

			if(sp.size() <= 2) continue;

			vector<int> vv(sp.begin(), sp.end());
			hyper_edge he(vv, 1);
			map<hyper_edge, int>::iterator it = mhe.find(he);
			if(it == mhe.end()) mhe.insert(pair<hyper_edge, int>(he, 1));
			else it->second++;
		}
	}

	hedges.clear();
	for(map<hyper_edge, int>::iterator it = mhe.begin(); it != mhe.end(); it++)
	{
		hyper_edge he = it->first;
		he.count = it->second;
		hedges.push_back(he);
	}

	return 0;
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
		b.lrgn = li->second;
		b.rrgn = ri->second;
	}
	return 0;
}

int bundle::size() const
{
	return pexons.size();
}

int bundle::build_splice_graph(splice_graph &gr) const
{
	//if(pexons.size() == 0) return 0;
	// vertices: start, each region, end
	gr.add_vertex();
	gr.set_vertex_string(0, "");
	gr.set_vertex_weight(0, 0);
	gr.set_vertex_info(0, vertex_info());
	for(int i = 0; i < pexons.size(); i++)
	{
		const partial_exon &r = pexons[i];
		int length = r.rpos - r.lpos;
		assert(length >= 1);
		gr.add_vertex();
		gr.set_vertex_string(i + 1, r.label());
		gr.set_vertex_weight(i + 1, r.ave);
		vertex_info vi;
		vi.length = length;
		vi.stddev = r.dev < 1.0 ? 1.0 : r.dev;
		vi.adjust = r.adjust;
		gr.set_vertex_info(i + 1, vi);
	}

	gr.add_vertex();
	gr.set_vertex_string(pexons.size() + 1, "");
	gr.set_vertex_weight(pexons.size() + 1, 0);
	gr.set_vertex_info(pexons.size() + 1, vertex_info());

	// edges: connecting adjacent pexons => e2w
	for(int i = 0; i < (int)(pexons.size()) - 1; i++)
	{
		const partial_exon &x = pexons[i];
		const partial_exon &y = pexons[i + 1];

		if(x.rpos != y.lpos) continue;

		assert(x.rpos == y.lpos);
		int32_t xr = compute_overlap(imap, x.rpos - 1);
		int32_t yl = compute_overlap(imap, y.lpos);
		// TODO
		double wt = xr < yl ? xr : yl;
		double sd = 1.0;
		//double sd = 0.5 * x.dev + 0.5 * y.dev;

		edge_descriptor p = gr.add_edge(i + 1, i + 2);
		gr.set_edge_weight(p, wt);
		gr.set_edge_info(p, edge_info());
	}

	// edges: each junction => and e2w
	for(int i = 0; i < junctions.size(); i++)
	{
		const junction &b = junctions[i];
		const partial_exon &x = pexons[b.lrgn];
		const partial_exon &y = pexons[b.rrgn];

		//double sd = 0.5 * x.dev + 0.5 * y.dev;
		double sd = 1.0;
		edge_descriptor p = gr.add_edge(b.lrgn + 1, b.rrgn + 1);
		gr.set_edge_weight(p, b.count);
		gr.set_edge_info(p, edge_info());
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
			gr.set_edge_weight(p, 0);
			gr.set_edge_info(p, edge_info());
		}

		if(r.rtype == END_BOUNDARY) 
		{
			edge_descriptor p = gr.add_edge(i + 1, tt);
			gr.set_edge_weight(p, 0);
			gr.set_edge_info(p, edge_info());
		}
	}

	return 0;
}

int bundle::print(int index) const
{
	printf("\nBundle %d: ", index);
	printf("tid = %d, #hits = %lu, range = %s:%d-%d, orient = %c %s\n",
			tid, hits.size(), chrm.c_str(), lpos, rpos, strand, (pexons.size() >= 100 ? "[monster]" : ""));
	// print hits
	/*
	for(int i = 0; i < hits.size(); i++)
	{
		hits[i].print();
	}
	*/

	// print super regions
	/*
	for(int i = 0; i < srs.size(); i++)
	{
		srs[i].print(i);
	}
	*/

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

	for(int i = 0; i < hedges.size(); i++)
	{
		hedges[i].print(i);
	}

	return 0;
}

int bundle::output_gtf(ofstream &fout, const vector<path> &paths, const string &prefix, int index) const
{
	fout.precision(2);
	fout<<fixed;

	for(int i = 0; i < paths.size(); i++)
	{
		const vector<int> &v = paths[i].v;
		double abd = paths[i].abd;
		if(v.size() < 2) continue;

		if(abd <= transcript_min_expression) continue;

		fout<<chrm.c_str()<<"\t";		// chromosome name
		fout<<prefix.c_str()<<"\t";		// source
		fout<<"transcript\t";			// feature
		fout<<lpos<<"\t";				// left position
		fout<<rpos<<"\t";				// right position
		fout<<1000<<"\t";				// score, now as abundance
		fout<<strand<<"\t";				// strand
		fout<<".\t";					// frame
		fout<<"gene_id \""<<algo.c_str()<<"."<<index<<"\"; ";
		fout<<"transcript_id \""<<algo.c_str()<<"."<<index<<"."<<i + 1<<"\"; ";
		fout<<"expression \""<<abd<<"\";"<<endl;

		assert(v[0] == 0);
		assert(v[v.size() - 1] == pexons.size() + 1);

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
			fout<<prefix.c_str()<<"\t";			// source
			fout<<"exon\t";						// feature
			fout<<lower(it->first) + 1<<"\t";	// left position
			fout<<upper(it->first)<<"\t";		// right position
			fout<<1000<<"\t";					// score
			fout<<strand<<"\t";					// strand
			fout<<".\t";						// frame
			fout<<"gene_id \""<<algo.c_str()<<"."<<index<<"\"; ";
			fout<<"transcript_id \""<<algo.c_str()<<"."<<index<<"."<<i + 1<<"\"; ";
			fout<<"exon_number \""<<cnt<<"\"; ";
			fout<<"expression \""<<abd<<"\";"<<endl;
		}
	}
	return 0;
}

int bundle::print_5end_coverage() const
{
	if(pexons.size() == 0) return 0;

	int n = 500;
	const partial_exon &p = pexons[0];
	if(p.rpos - p.lpos <= n) return 0;

	int k = 0;
	for(int32_t i = p.lpos; i < p.lpos + n; i++)
	{
		int32_t c = compute_overlap(imap, i);
		printf("5end coverage %d %d\n", k++, c);
	}

	return 0;
}

int bundle::print_3end_coverage() const
{
	if(pexons.size() == 0) return 0;

	int n = 500;
	const partial_exon &p = pexons[pexons.size() - 1];
	if(p.rpos - p.lpos <= n) return 0;

	int k = 0;
	for(int32_t i = p.rpos - 1; i >= p.rpos - n; i--)
	{
		int32_t c = compute_overlap(imap, i);
		printf("3end coverage %d %d\n", k++, c);
	}

	return 0;
}
