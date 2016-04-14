#include <cassert>
#include <cstdio>
#include <map>
#include <iomanip>
#include <fstream>

#include "bundle.h"
#include "binomial.h"

bundle::bundle(const bundle_base &bb)
	: bundle_base(bb)
{
}

bundle::~bundle()
{}

int bundle::build()
{
	// make sure all reads are sorted 
	check_left_ascending();

	build_split_interval_map();

	infer_junctions();

	build_partial_exons();

	link_regions();
	split_region_boundaries();

	return 0;
}

int bundle::build_split_interval_map()
{
	imap.clear();
	vector<int64_t> v;
	for(int i = 0; i < hits.size(); i++)
	{
		hits[i].get_matched_intervals(v);
		for(int k = 0; k < v.size(); k++)
		{
			int32_t s = high32(v[k]);
			int32_t t = low32(v[k]);
			imap += make_pair(ROI(s, t), 1);
		}
	}
	return 0;
}

int bundle::locate_hits(int32_t p, int &li)
{
	li = -1;
	hit h(p);
	vector<hit>::iterator low = lower_bound(hits.begin(), hits.end(), h);
	vector<hit>::iterator up = upper_bound(hits.begin(), hits.end(), h);
	if(low == hits.end()) return 0;
	li = low - hits.begin();
	return (up - low);
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

int bundle::add_start_boundary()
{
	for(int i = 0; i < boundaries.size(); i++)
	{
		boundary &b = boundaries[i];
		if(b.pos == lpos) return 0;
	}

	boundary b(START_BOUNDARY, lpos, 1, 100, 100);
	boundaries.push_back(b);

	return 0;
}

int bundle::add_end_boundary()
{
	for(int i = 0; i < boundaries.size(); i++)
	{
		boundary &b = boundaries[i];
		if(b.pos == rpos) return 0;
	}

	boundary b(END_BOUNDARY, rpos, 1, 100, 100);
	boundaries.push_back(b);

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

int bundle::build_regions()
{
	typedef map<int32_t, int> MM;
	typedef pair<int32_t, int> PP;

	MM s;
	for(int i = 0; i < junctions.size(); i++)
	{
		s.insert(PP(junctions[i].lpos, LEFT_SPLICE));
		s.insert(PP(junctions[i].rpos, RIGHT_SPLICE));
	}

	for(int i = 0; i < boundaries.size(); i++)
	{
		s.insert(PP(boundaries[i].pos, boundaries[i].type));
	}

	vector<int32_t> v;
	MM::iterator it;
	for(it = s.begin(); it != s.end(); it++) v.push_back(it->first);
	sort(v.begin(), v.end());

	if(v.size() <= 1) return 0;

	for(int i = 0; i < v.size() - 1; i++)
	{
		region r(v[i], v[i + 1], s[v[i]], s[v[i + 1]], &imap);
		regions.push_back(r);
	}

	for(int i = 0; i < regions.size(); i++)
	{
		regions[i].build();
	}
	return 0;
}

int bundle::link_regions()
{
	MPI lm;
	MPI rm;
	for(int i = 0; i < regions.size(); i++)
	{
		int32_t l = regions[i].lpos;
		int32_t r = regions[i].rpos;
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

int bundle::split_region_boundaries()
{
	for(int i = 0; i < regions.size(); i++)
	{
		create_split(imap, regions[i].lpos);
		create_split(imap, regions[i].rpos);
	}
	return 0;
}

int bundle::print(int index) const
{
	printf("\nBundle %d: ", index);
	char o = phits > qhits ? '+' : '-';
	printf("tid = %d, #hits = %lu (%d, %d), range = %s:%d-%d, orient = %c\n",
			tid, hits.size(), phits, qhits, chrm.c_str(), lpos, rpos, o);
	// print hits
	/*
	for(int i = 0; i < hits.size(); i++)
	{
		hits[i].print();
	}
	*/

	// print junctions 
	for(int i = 0; i < junctions.size(); i++)
	{
		junctions[i].print(i);
	}

	// print boundaries
	for(int i = 0; i < boundaries.size(); i++)
	{
		boundaries[i].print(i);
	}

	// print regions
	for(int i = 0; i < regions.size(); i++)
	{
		regions[i].print(i);
	}

	// print region boundaries
	for(int i = 0; i < regions.size(); i++)
	{
		regions[i].print_boundaries(i);
	}

	return 0;
}

int bundle::build_splice_graph(splice_graph &gr) const
{
	// vertices: start, each region, end
	gr.add_vertex();
	gr.set_vertex_string(0, "");
	gr.set_vertex_weight(0, 0);
	gr.set_vertex_stddev(0, 1);
	for(int i = 0; i < regions.size(); i++)
	{
		const region &r = regions[i];
		double ave = 0;
		double dev = 1;
		if(r.empty == false)
		{
			ave = r.ave_abd;
			dev = r.dev_abd;					
			//dev = r.dev_abd / sqrt(r.rcore - r.lcore); // TODO
		}
		gr.add_vertex();
		gr.set_vertex_string(i + 1, r.label());
		gr.set_vertex_weight(i + 1, ave);
		gr.set_vertex_stddev(i + 1, dev);
	}
	gr.add_vertex();
	gr.set_vertex_string(regions.size() + 1, "");
	gr.set_vertex_weight(regions.size() + 1, 0);
	gr.set_vertex_stddev(regions.size() + 1, 1);

	// edges: connecting adjacent regions => e2w
	for(int i = 0; i < regions.size() - 1; i++)
	{
		const region &x = regions[i];
		const region &y = regions[i + 1];

		if(x.empty || y.empty) continue;

		if(x.right_break()) continue;
		if(y.left_break()) continue;

		//if(x.rtype == RIGHT_BOUNDARY) continue;
		//if(y.ltype == LEFT_BOUNDARY) continue;

		assert(x.rpos == y.lpos);
		int32_t xr = compute_overlap(imap, x.rpos - 1);
		int32_t yl = compute_overlap(imap, y.lpos);
		double wt = xr < yl ? xr : yl;
		double sd = 0.5 * x.dev_abd + 0.5 * y.dev_abd;

		edge_descriptor p = gr.add_edge(i + 1, i + 2);
		gr.set_edge_weight(p, wt);
		gr.set_edge_stddev(p, sd);
	}

	// edges: each junction => and e2w
	for(int i = 0; i < junctions.size(); i++)
	{
		const junction &b = junctions[i];
		const region &x = regions[b.lrgn];
		const region &y = regions[b.rrgn];

		double sd = 0.5 * x.dev_abd + 0.5 * y.dev_abd;

		edge_descriptor p = gr.add_edge(b.lrgn + 1, b.rrgn + 1);
		gr.set_edge_weight(p, b.count);
		gr.set_edge_stddev(p, sd);
	}

	// edges: connecting start/end and regions
	int ss = 0;
	int tt = regions.size() + 1;
	for(int i = 0; i < regions.size(); i++)
	{
		const region &r = regions[i];
		if(r.empty == true) continue;

		// TODO
		//if(r.ltype == LEFT_BOUNDARY || r.ltype == START_BOUNDARY)
		if(r.left_break() || r.ltype == LEFT_BOUNDARY || r.ltype == START_BOUNDARY)
		{
			edge_descriptor p = gr.add_edge(ss, i + 1);
			gr.set_edge_weight(p, r.ave_abd);
			gr.set_edge_stddev(p, r.dev_abd);
		}

		// TODO
		//if(r.rtype == RIGHT_BOUNDARY || r.rtype == END_BOUNDARY) 
		if(r.right_break() || r.rtype == RIGHT_BOUNDARY || r.rtype == END_BOUNDARY) 
		{
			edge_descriptor p = gr.add_edge(i + 1, tt);
			gr.set_edge_weight(p, r.ave_abd);
			gr.set_edge_stddev(p, r.dev_abd);
		}
	}

	// check
	if(regions.size() == 1) return 0;
	for(int i = 0; i < regions.size(); i++)
	{
		const region &r = regions[i];
		if(r.empty == false) continue;
		assert(gr.in_degree(i + 1) == 0);
		assert(gr.out_degree(i + 1) == 0);
		// TODO, still bugy here
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
		fout<<lpos<<"\t";				// left position
		fout<<rpos<<"\t";				// right position
		fout<<1000<<"\t";				// score, now as abundance
		fout<<"+\t";					// strand
		fout<<".\t";					// frame
		fout<<"gene_id \""<<algo.c_str()<<"."<<index<<"\"; ";
		fout<<"transcript_id \""<<algo.c_str()<<"."<<index<<"."<<i + 1<<"\"; ";
		fout<<"expression \""<<abd<<"\";"<<endl;

		assert(v[0] == 0);
		assert(v[v.size() - 1] == regions.size() + 1);

		join_interval_map jmap;
		for(int k = 1; k < v.size() - 1; k++)
		{
			const region &r = regions[v[k] - 1];
			jmap += make_pair(ROI(r.lpos, r.rpos), 1);
		}

		int cnt = 0;
		for(JIMI it = jmap.begin(); it != jmap.end(); it++)
		{
			fout<<chrm.c_str()<<"\t";			// chromosome name
			fout<<prefix.c_str()<<"\t";			// source
			fout<<"exon\t";						// feature
			fout<<lower(it->first) + 1<<"\t";		// left position
			fout<<upper(it->first)<<"\t";	// right position
			fout<<1000<<"\t";					// score
			fout<<"+\t";						// strand
			fout<<".\t";						// frame
			fout<<"gene_id \""<<algo.c_str()<<"."<<index<<"\"; ";
			fout<<"transcript_id \""<<algo.c_str()<<"."<<index<<"."<<i + 1<<"\"; ";
			fout<<"exon_number \""<<cnt<<"\"; ";
			fout<<"expression \""<<abd<<"\";"<<endl;
		}
	}
	return 0;
}


