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
	check_left_ascending();
	build_split_interval_map();
	infer_junctions();
	build_partial_exons();
	link_partial_exons();
	//split_boundaries();

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


int bundle::build_partial_exons()
{
	vector<PPI> v;
	v.push_back(PPI(lpos, START_BOUNDARY));
	v.push_back(PPI(rpos, END_BOUNDARY));

	for(int i = 0; i < junctions.size(); i++)
	{
		v.push_back(PPI(junctions[i].lpos, LEFT_SPLICE));
		v.push_back(PPI(junctions[i].rpos, RIGHT_SPLICE));
	}

	sort(v.begin(), v.end());

	for(int i = 0; i < v.size() - 1; i++)
	{
		//printf("try region %d\n", i);
		region r(v[i].first, v[i + 1].first, v[i].second, v[i + 1].second, &imap);
		vector<partial_exon> vp = r.build();
		//r.print(i);
		pexons.insert(pexons.end(), vp.begin(), vp.end());
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

/*
int bundle::split_boundaries()
{
	for(int i = 0; i < pexons.size(); i++)
	{
		create_split(imap, pexons[i].lpos);
		create_split(imap, pexons[i].rpos);
	}
	return 0;
}
*/

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
		gr.set_vertex_weight(i + 1, r.ave_abd);
		gr.set_vertex_info(i + 1, vertex_info(length));
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
		double wt = xr < yl ? xr : yl;
		double sd = 1.0;
		//double sd = 0.5 * x.dev_abd + 0.5 * y.dev_abd;

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

		//double sd = 0.5 * x.dev_abd + 0.5 * y.dev_abd;
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
	char o = phits > qhits ? '+' : '-';
	printf("tid = %d, #hits = %lu (%d, %d), range = %s:%d-%d, orient = %c %s\n",
			tid, hits.size(), phits, qhits, chrm.c_str(), lpos, rpos, o, (pexons.size() >= 100 ? "[monster]" : ""));
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

	// print partial exons
	for(int i = 0; i < pexons.size(); i++)
	{
		pexons[i].print(i);
	}

	return 0;
}

int bundle::output_gtf(ofstream &fout, const vector<path> &paths, const string &prefix, int index) const
{
	fout.precision(2);
	fout<<fixed;

	char o = phits > qhits ? '+' : '-';

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
		fout<<o<<"\t";					// strand
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
			fout<<o<<"\t";						// strand
			fout<<".\t";						// frame
			fout<<"gene_id \""<<algo.c_str()<<"."<<index<<"\"; ";
			fout<<"transcript_id \""<<algo.c_str()<<"."<<index<<"."<<i + 1<<"\"; ";
			fout<<"exon_number \""<<cnt<<"\"; ";
			fout<<"expression \""<<abd<<"\";"<<endl;
		}
	}
	return 0;
}
