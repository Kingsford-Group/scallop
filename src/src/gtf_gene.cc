#include <algorithm>
#include "gtf_gene.h"

int gtf_gene::build_splice_graph(splice_graph &gr)
{
	sort(exons.begin(), exons.end());
	build_transcripts();
	build_split_interval_map();
	add_vertices(gr);
	add_edges(gr);
	return 0;
}

int gtf_gene::add_vertices(splice_graph &gr)
{
	gr.add_vertex();
	int32_t s = compute_sum_expression();
	gr.set_vertex_weight(0, s);
	gr.set_vertex_stddev(0, 1.0);

	SIMI it;
	for(it = imap.begin(); it != imap.end(); it++)
	{
		gr.add_vertex();
		gr.set_vertex_weight(gr.num_vertices() - 1, it->second);
		gr.set_vertex_stddev(gr.num_vertices() - 1, 1.0);
	}

	gr.add_vertex();
	gr.set_vertex_weight(gr.num_vertices() - 1, s);
	gr.set_vertex_stddev(gr.num_vertices() - 1, 1.0);
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
			SIMI it = imap.find(ge.start);
			assert(it != imap.end());
			while(true)
			{
				int uu = distance((SIMI)(imap.begin()), it) + 1;
				add_single_edge(u, uu, expr, gr);
				u = uu;
				if(upper(it->first) >= ge.end) break;
				it++;
			}
		}
		add_single_edge(u, gr.num_vertices() -1, expr, gr);
	}
	return 0;
}

int gtf_gene::add_single_edge(int s, int t, double w, splice_graph &gr)
{
	PEB p = gr.edge(s, t);
	if(p.second == true)
	{
		double w0 = gr.get_edge_weight(p.first);	
		gr.set_edge_weight(p.first, w + w0);
	}
	else
	{
		edge_descriptor p = gr.add_edge(s, t);
		gr.set_edge_weight(p, w);
		gr.set_edge_stddev(p, 1.0);
	}
	return 0;
}

int gtf_gene::build_split_interval_map()
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

	for(SIMI it = imap.begin(); it != imap.end(); it++)
	{
		printf("%lu vertex: [%d, %d)\n", distance((SIMI)(imap.begin()), it), lower(it->first), upper(it->first));
	}
	return 0;
}


int gtf_gene::output_gtf(ofstream &fout, const vector<path> &paths, const string &prefix) const
{
	fout.precision(2);
	fout<<fixed;

	if(exons.size() == 0) return 0;

	string chrm = exons[0].seqname;
	string gene = exons[0].gene_id;

	for(int i = 0; i < paths.size(); i++)
	{
		const vector<int> &v = paths[i].v;
		double abd = paths[i].abd;
		if(v.size() < 2) continue;

		SIMI si = imap.begin();
		SIMI ti = imap.end();
		ti--;

		fout<<chrm.c_str()<<"\t";					// chromosome name
		fout<<prefix.c_str()<<"\t";					// source
		fout<<"transcript\t";						// feature
		fout<<lower(si->first)<<"\t";				// left position
		fout<<upper(ti->first)<<"\t";				// right position
		fout<<1000<<"\t";							// score, now as abundance
		fout<<"+\t";								// strand
		fout<<".\t";								// frame
		fout<<"gene_id \""<<gene.c_str()<<"\"; ";
		fout<<"transcript_id \""<<gene.c_str()<<"."<<i + 1<<"\"; ";
		fout<<"abundance \""<<abd<<"\";"<<endl;

		assert(v[0] == 0);
		join_interval_map jmap;
		for(int k = 1; k < v.size() - 1; k++)
		{
			SIMI it = imap.begin();
			advance(it, v[k] - 1);
			jmap += make_pair(ROI(lower(it->first), upper(it->first)), 1);
		}

		int cnt = 0;
		for(JIMI it = jmap.begin(); it != jmap.end(); it++)
		{
			fout<<chrm.c_str()<<"\t";			// chromosome name
			fout<<prefix.c_str()<<"\t";			// source
			fout<<"exon\t";						// feature
			fout<<lower(it->first)<<"\t";		// left position
			fout<<upper(it->first) - 1<<"\t";	// right position
			fout<<1000<<"\t";					// score, now as abundance
			fout<<"+\t";						// strand
			fout<<".\t";						// frame
			fout<<"gene_id \""<<gene.c_str()<<"\"; ";
			fout<<"transcript_id \""<<gene.c_str()<<"."<<i + 1<<"\"; ";
			fout<<"exon \""<<++cnt<<"\";"<<endl;
		}
	}
	return 0;
}

int gtf_gene::output_gtf(ofstream &fout) const
{
	fout.precision(2);
	fout<<fixed;

	if(exons.size() == 0) return 0;

	string chrm = exons[0].seqname;
	string gene = exons[0].gene_id;
	string src = exons[0].source;

	for(int i = 0; i < transcripts.size(); i++)
	{
		const vector<int> &v = transcripts[i];
		if(v.size() == 0) continue;

		double abd = exons[v[0]].expression;

		int lpos = exons[v[0]].start;
		int rpos = exons[v[v.size() - 1]].end;

		fout<<chrm.c_str()<<"\t";					// chromosome name
		fout<<src.c_str()<<"\t";					// source
		fout<<"transcript\t";						// feature
		fout<<lpos<<"\t";							// left position
		fout<<rpos<<"\t";							// right position
		fout<<1000<<"\t";							// score, now as abundance
		fout<<"+\t";								// strand
		fout<<".\t";								// frame
		fout<<"gene_id \""<<gene.c_str()<<"\"; ";
		fout<<"transcript_id \""<<gene.c_str()<<"."<<i + 1<<"\"; ";
		fout<<"abundance \""<<abd<<"\";"<<endl;

		for(int k = 0; k < v.size(); k++)
		{
			const gtf_exon &ge = exons[v[k]];
			fout<<chrm.c_str()<<"\t";			// chromosome name
			fout<<src.c_str()<<"\t";			// source
			fout<<"exon\t";						// feature
			fout<<ge.start<<"\t";				// left position
			fout<<ge.end - 1<<"\t";					// right position
			fout<<1000<<"\t";					// score, now as abundance
			fout<<"+\t";						// strand
			fout<<".\t";						// frame
			fout<<"gene_id \""<<gene.c_str()<<"\"; ";
			fout<<"transcript_id \""<<gene.c_str()<<"."<<i + 1<<"\"; ";
			fout<<"exon \""<<k + 1<<"\";"<<endl;
		}
	}
	return 0;
}
