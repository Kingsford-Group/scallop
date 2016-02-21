#include <algorithm>
#include "gtf_gene.h"

int gtf_gene::build_splice_graph(splice_graph &gr)
{
	sort(exons.begin(), exons.end());
	build_transcripts();
	build_interval_map();
	add_vertices(gr);
	add_edges(gr);
	return 0;
}

int gtf_gene::add_vertices(splice_graph &gr)
{
	add_vertex(gr);
	int32_t s = compute_sum_expression();
	put(get(vertex_weight, gr), 0, s);
	put(get(vertex_stddev, gr), 0, 1.0);

	ICI it;
	for(it = imap.begin(); it != imap.end(); it++)
	{
		add_vertex(gr);
		put(get(vertex_weight, gr), num_vertices(gr) - 1, it->second);
		put(get(vertex_stddev, gr), num_vertices(gr) - 1, 1.0);
	}

	add_vertex(gr);
	put(get(vertex_weight, gr), num_vertices(gr) - 1, s);
	put(get(vertex_stddev, gr), num_vertices(gr) - 1, 1.0);
	return 0;
}

int gtf_gene::add_edges(splice_graph &gr)
{
	for(int i = 0; i < transcripts.size(); i++)
	{
		vector<int> &v = transcripts[i];
		assert(v.size() >= 1);
		int32_t expr = exons[v[0]].expression;
		int u = 0;
		for(int k = 0; k < v.size(); k++)
		{
			gtf_exon &ge = exons[v[k]];
			ICI it = imap.find(ge.start);
			assert(it != imap.end());
			while(true)
			{
				int uu = distance((ICI)(imap.begin()), it) + 1;
				add_single_edge(u, uu, expr, gr);
				u = uu;
				if(upper(it->first) >= ge.end) break;
				it++;
			}
		}
		add_single_edge(u, num_vertices(gr) -1, expr, gr);
	}
	return 0;
}

int gtf_gene::add_single_edge(int s, int t, double w, splice_graph &gr)
{
	PEB p = edge(s, t, gr);
	if(p.second == true)
	{
		double w0 = get(get(edge_weight, gr), p.first);	
		put(get(edge_weight, gr), p.first, w + w0);
	}
	else
	{
		PEB p = add_edge(s, t, gr);
		assert(p.second == true);
		put(get(edge_weight, gr), p.first, w);
		put(get(edge_stddev, gr), p.first, 1.0);
	}
	return 0;
}

int gtf_gene::build_interval_map()
{
	imap.clear();
	for(int i = 0; i < exons.size(); i++)
	{
		gtf_exon &ge = exons[i];
		imap += make_pair(ROI(ge.start, ge.end), ge.expression);
	}
	return 0;
}

int gtf_gene::add_exon(const gtf_exon &ge)
{
	exons.push_back(ge);
	return 0;
}

int gtf_gene::build_transcripts()
{
	map<string, int> m;
	for(int i = 0; i < exons.size(); i++)
	{
		gtf_exon &ge = exons[i];
		if(m.find(ge.transcript_id) == m.end())
		{
			vector<int> v;
			v.push_back(i);
			m.insert(pair<string, int>(ge.transcript_id, transcripts.size()));
			transcripts.push_back(v);
		}
		else
		{
			transcripts[m[ge.transcript_id]].push_back(i);
		}
	}
	return 0;
}

int32_t gtf_gene::compute_sum_expression()
{
	int32_t s = 0;
	for(int i = 0; i < transcripts.size(); i++)
	{
		vector<int> & v = transcripts[i];
		assert(v.size() >= 1);
		s += exons[v[0]].expression;
	}
	return s;
}

int gtf_gene::print()
{
	for(int i = 0; i < exons.size(); i++)
	{
		exons[i].print();
	}
	return 0;
}
