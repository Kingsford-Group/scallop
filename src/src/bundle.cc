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
	build_junctions();
	build_junction_graph();
	align_hits();
	build_regions();
	build_partial_exons();
	build_hyper_edges();
	link_partial_exons();
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
		for(int k = 0; k < v.size(); k++)
		{
			int64_t p = v[k];
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

int bundle::build_junction_graph()
{
	MPI s;
	s.insert(PI(lpos, START_BOUNDARY));
	s.insert(PI(rpos, END_BOUNDARY));
	for(int i = 0; i < junctions.size(); i++)
	{
		int32_t l = junctions[i].lpos;
		int32_t r = junctions[i].rpos;

		if(s.find(l) == s.end()) s.insert(PI(l, LEFT_SPLICE));
		else if(s[l] == RIGHT_SPLICE) s[l] = LEFT_RIGHT_SPLICE;

		if(s.find(r) == s.end()) s.insert(PI(r, RIGHT_SPLICE));
		else if(s[r] == LEFT_SPLICE) s[r] = LEFT_RIGHT_SPLICE;
	}

	vector<PPI> v(s.begin(), s.end());
	sort(v.begin(), v.end());

	// position to vertex
	MPI p2v;
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
		vi.pos = v[i].first;
		vi.type = v[i].second;
		jr.set_vertex_info(i, vi);
	}

	// edges: connecting adjacent regions
	for(int i = 0; i < v.size() - 1; i++)
	{
		edge_descriptor p = jr.add_edge(i, i + 1);
		edge_info ei;
		ei.jid = -1;
		jr.set_edge_info(p, ei);
	}

	// edges: each junction
	for(int i = 0; i < junctions.size(); i++)
	{
		int32_t l = junctions[i].lpos;
		int32_t r = junctions[i].rpos;
		assert(p2v.find(l) != p2v.end());
		assert(p2v.find(r) != p2v.end());
		edge_descriptor p = jr.add_edge(p2v[l], p2v[r]);
		edge_info ei;
		ei.jid = i;
		jr.set_edge_info(p, ei);
	}

	return 0;
}

int bundle::draw_junction_graph(const string &file)
{
	char buf[10240];

	MIS mis;
	for(int i = 0; i < jr.num_vertices(); i++)
	{
		int32_t p = jr.get_vertex_info(i).pos;
		sprintf(buf, "%d", p);
		mis.insert(PIS(i, buf));
	}

	MES mes;
	edge_iterator it1, it2;
	for(tie(it1, it2) = jr.edges(); it1 != it2; it1++)
	{
		int jid = jr.get_edge_info(*it1).jid;
		sprintf(buf, "%d", jid);
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
				int id = jr.get_edge_info(ve[k]).jid;
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
		int32_t p1 = jr.get_vertex_info(m).pos;
		int32_t p2 = jr.get_vertex_info(m + 1).pos;
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
			int id = jr.get_edge_info(*it1).jid;
			if(id < 0) 
			{
				assert(ss == k - 1);
				w = jr.get_vertex_info(k).pos - jr.get_vertex_info(ss).pos;

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
		int id = jr.get_edge_info(ve2[k]).jid;
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
			int id = jr.get_edge_info(*it1).jid;
			if(id < 0) 
			{
				assert(ss == k - 1);
				w = jr.get_vertex_info(k).pos - jr.get_vertex_info(ss).pos;
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

int bundle::align_hits()
{
	imap.clear();
	jmap.clear();
	for(int i = 0; i < hits.size(); i++)
	{
		hit &h = hits[i];
		vector<int64_t> vm;
		vector<int64_t> vi;
		vector<int64_t> vd;

		h.get_mid_intervals(vm, vi, vd);

		h.print();

		for(int k = 0; k < vm.size(); k++)
		{
			int32_t s = high32(vm[k]);
			int32_t t = low32(vm[k]);
			imap += make_pair(ROI(s, t), 1);
		}

		for(int k = 0; k < vi.size(); k++)
		{
			int32_t s = high32(vi[k]);
			int32_t t = low32(vi[k]);
			printf("add insertion (%d, %d)\n", s, t);
			jmap += make_pair(ROI(s, t), 1);
		}

		for(int k = 0; k < vd.size(); k++)
		{
			int32_t s = high32(vd[k]);
			int32_t t = low32(vd[k]);
			printf("add deletion (%d, %d)\n", s, t);
			jmap += make_pair(ROI(s, t), 1);
		}
	}

	// Test jmap
	SIMI::iterator jit;
	for(jit = jmap.begin(); jit != jmap.end(); jit++)
	{
		printf("jmap (%d, %d) -> %d\n", lower(jit->first), upper(jit->first), jit->second);
	}

	return 0;
}

bool bundle::verify_unique_mapping(const hit &h)
{
	if((h.flag & 0x8) >= 1) return true;
	if(h.mpos < h.rpos) return true;

	int li = search_junction_graph(h.rpos - 1);
	if(li < 0) return false;
	assert(li >= 0 && li < jr.num_vertices() - 1);

	int32_t l1 = jr.get_vertex_info(li).pos;
	int32_t l2 = jr.get_vertex_info(li + 1).pos;
	assert(h.rpos - 1 >= l1 && h.rpos - 1 < l2);
	
	//printf("li = %d, l1 = %d, l2 = %d\n", li, l1, l2);
	if(h.mpos < l2) return true;

	int ri = search_junction_graph(h.mpos);
	if(ri < 0) return false;

	assert(ri >= 0 && ri < jr.num_vertices() - 1);
	assert(ri > li);

	int32_t r1 = jr.get_vertex_info(ri).pos;
	int32_t r2 = jr.get_vertex_info(ri + 1).pos;
	assert(h.mpos >= r1 && h.mpos < r2);

	int d1 = r1 - l2;
	assert(d1 >= 0);

	int d2 = traverse_junction_graph1(li + 1, ri);

	//printf("ri = %d, r1 = %d, r2 = %d\n", ri, r1, r2);
	//printf("d1 = %d, d2 = %d\n", d1, d2);

	if(d2 <= -1 || d1 < d2) return true;
	else return false;
}

int bundle::compute_read1_intervals(const hit &h, vector<int64_t> &vv)
{
	assert(verify_unique_mapping(h) == true);

	vv.clear();
	if(h.isize < 0) return 0;

	vector<int64_t> vm;
	h.get_matched_intervals(vm);
	for(int k = 0; k < vm.size(); k++)
	{
		int32_t s = high32(vm[k]);
		int32_t t = low32(vm[k]);
		if(s >= h.mpos) break;
		if(t >= h.mpos) t = h.mpos;
		vv.push_back(vm[k]);
	}

	if(h.rpos >= h.mpos) return 0;

	int li = search_junction_graph(h.rpos - 1);
	assert(li >= 0 && li < jr.num_vertices() - 1);
	int32_t l1 = jr.get_vertex_info(li).pos;
	int32_t l2 = jr.get_vertex_info(li + 1).pos;
	assert(h.rpos - 1 >= l1 && h.rpos - 1 < l2);
	
	if(h.mpos <= l2)
	{
		vv.push_back(pack(h.rpos, h.mpos));
		return 0;
	}

	int ri = search_junction_graph(h.mpos);
	assert(ri >= 0 && ri < jr.num_vertices() - 1);
	assert(ri > li);
	int32_t r1 = jr.get_vertex_info(ri).pos;
	int32_t r2 = jr.get_vertex_info(ri + 1).pos;
	assert(h.mpos >= r1 && h.mpos < r2);

	vv.push_back(pack(h.rpos, l2));
	vv.push_back(pack(r1, h.mpos));
	return 0;
}

int bundle::build_regions()
{
	for(int k = 0; k < jr.num_vertices() - 1; k++)
	{
		int32_t lpos = jr.get_vertex_info(k).pos;
		int32_t rpos = jr.get_vertex_info(k + 1).pos;

		int ltype = jr.get_vertex_info(k).type; 
		int rtype = jr.get_vertex_info(k + 1).type; 

		if(ltype == LEFT_RIGHT_SPLICE) ltype = RIGHT_SPLICE;
		if(rtype == LEFT_RIGHT_SPLICE) rtype = LEFT_SPLICE;

		regions.push_back(region(lpos, rpos, ltype, rtype, &imap, &jmap));
	}

	return 0;
}

int bundle::build_partial_exons()
{
	pexons.clear();
	for(int i = 0; i < regions.size(); i++)
	{
		region &r = regions[i];
		vector<partial_exon> v;
		r.build_partial_exons(v);
		pexons.insert(pexons.end(), v.begin(), v.end());
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

int bundle::build_hyper_edges()
{
	vector<int> list;
	map<string, PI> m;
	vector<int64_t> v;
	map<hyper_edge, int> mhe;
	for(int i = 0; i < hits.size(); i++)
	{
		hit &h = hits[i];

		if(h.isize >= 0 && verify_unique_mapping(h) == false) continue;

		if(h.isize < 0) h.get_matched_intervals(v);
		else compute_read1_intervals(h, v);

		if(v.size() == 0) continue;
		string s = h.qname;

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
			}

			/*
			h.print();
			printf(" -> ");
			for(int k = pi.first; k < pi.second; k++) printf(" %d", list[k]);
			printf("\n");
			*/

			if(pi.second > pi.first)
			{
				m.insert(pair<string, PI>(s, pi));
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

			// hyper junctions

			if(sp.size() <= 2) continue;

			vector<int> vv(sp.begin(), sp.end());
			hyper_edge he(vv, 1);

			/*
			h.print();
			printf(" -> ");
			he.print(99);
			printf("\n");
			*/

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

int bundle::size() const
{
	return pexons.size();
}

int bundle::build_splice_graph(splice_graph &gr, vector<hyper_edge> &vhe) const
{
	//if(pexons.size() == 0) return 0;
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
		vi.stddev = r.dev < 1.0 ? 1.0 : r.dev;
		//vi.adjust = r.adjust;
		gr.set_vertex_info(i + 1, vi);
	}

	gr.add_vertex();
	vertex_info vin;
	vin.lpos = rpos;
	vin.rpos = rpos;
	gr.set_vertex_weight(pexons.size() + 1, 0);
	gr.set_vertex_info(pexons.size() + 1, vin);

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
		gr.set_edge_weight(p, wt < 1.0 ? 1.0 : wt);
		gr.set_edge_info(p, edge_info());
	}

	// edges: each junction => and e2w
	for(int i = 0; i < junctions.size(); i++)
	{
		const junction &b = junctions[i];

		if(b.lexon < 0 || b.rexon < 0) continue;

		const partial_exon &x = pexons[b.lexon];
		const partial_exon &y = pexons[b.rexon];

		//double sd = 0.5 * x.dev + 0.5 * y.dev;
		double sd = 1.0;
		edge_descriptor p = gr.add_edge(b.lexon + 1, b.rexon + 1);
		assert(b.count >= 1);
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

	// hyper-edges
	vhe.clear();
	for(int i = 0; i < hedges.size(); i++)
	{
		hyper_edge he = hedges[i];
		he.increase();
		vhe.push_back(he);
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

	// print hyper edges
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

		fout<<chrm.c_str()<<"\t";		// chromosome name
		fout<<prefix.c_str()<<"\t";		// source
		fout<<"transcript\t";			// feature
		fout<<lpos + 1<<"\t";			// left position
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
			fout<<"exon_number \""<<++cnt<<"\"; ";
			fout<<"expression \""<<abd<<"\";"<<endl;
		}
	}
	return 0;
}

