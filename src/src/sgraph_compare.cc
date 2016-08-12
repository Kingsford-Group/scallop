#include "sgraph_compare.h"
#include "util.h"
#include "draw.h"

sgraph_compare::sgraph_compare(const splice_graph &g1, const splice_graph &g2)
	:gr1(g1), gr2(g2)
{}

int sgraph_compare::compare(const string &file)
{
	imap.clear();
	build_split_interval_map(gr1);
	build_split_interval_map(gr2);

	add_vertices(gr3);
	add_inner_edges(gr1, gr3, 1);
	add_inner_edges(gr2, gr3, 2);
	add_existing_edges(gr1, gr3, 3);
	add_existing_edges(gr2, gr3, 4);

	draw(gr3, file);

	identify_unique_boundary_edges();

	return 0;
}

int sgraph_compare::build_split_interval_map(splice_graph &gr)
{
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		imap += make_pair(ROI(vi.lpos, vi.rpos), 1);
	}
	return 0;
}

int sgraph_compare::add_vertices(splice_graph &gr)
{
	if(imap.size() == 0) return 0;

	gr.add_vertex();
	gr.set_vertex_weight(0, 1);
	vertex_info vi0;
	SIMI it = imap.begin();
	vi0.lpos = lower(it->first);
	vi0.rpos = lower(it->first);
	gr.set_vertex_info(0, vi0);

	for(it = imap.begin(); it != imap.end(); it++)
	{
		gr.add_vertex();
		vertex_info vi;
		vi.length = upper(it->first) - lower(it->first);
		vi.lpos = lower(it->first);
		vi.rpos = upper(it->first);
		gr.set_vertex_weight(gr.num_vertices() - 1, 1);
		gr.set_vertex_info(gr.num_vertices() - 1, vi);
	}

	it = imap.end();
	it--;
	vertex_info vin;
	vin.lpos = upper(it->first);
	vin.rpos = upper(it->first);
	gr.add_vertex();
	gr.set_vertex_weight(gr.num_vertices() - 1, 1);
	gr.set_vertex_info(gr.num_vertices() - 1, vin);
	return 0;
}

int sgraph_compare::add_inner_edges(splice_graph &gt, splice_graph &gr, int type)
{
	for(int i = 1; i < gt.num_vertices() - 1; i++)
	{
		vertex_info vi = gt.get_vertex_info(i);
		double w = gt.get_vertex_weight(i);

		int k1 = search_splice_graph(gr, vi.lpos);
		int k2 = search_splice_graph(gr, vi.rpos - 1);
		assert(k1 >= 1 && k1 < gr.num_vertices() - 1);
		assert(k2 >= 1 && k2 < gr.num_vertices() - 1);

		for(int k = k1; k < k2; k++)
		{
			edge_descriptor p = gr.add_edge(k, k + 1);
			gr.set_edge_weight(p, w);
			edge_info ei;
			ei.type = type;
			gr.set_edge_info(p, ei);
		}
	}
	return 0;
}

int sgraph_compare::add_existing_edges(splice_graph &gt, splice_graph &gr, int type)
{
	edge_iterator it1, it2;
	for(tie(it1, it2) = gt.edges(); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		double w = gt.get_edge_weight(e);

		int s = e->source();
		int t = e->target();
		vertex_info vs = gt.get_vertex_info(s);
		vertex_info vt = gt.get_vertex_info(t);

		int ss = search_splice_graph(gr, vs.rpos - 1);
		int tt = search_splice_graph(gr, vt.lpos);

		//printf("edge = (%d, %d), type = %d, search %d -> %d\n", s, t, type, vs.rpos - 1, ss);
		//printf("edge = (%d, %d), type = %d, search %d -> %d\n", s, t, type, vt.lpos, tt);

		if(s == 0) ss = 0;
		if(t == gt.num_vertices() - 1) tt = gr.num_vertices() - 1;

		assert(ss >= 0 && ss < gr.num_vertices());
		assert(tt >= 0 && tt < gr.num_vertices());

		edge_descriptor p = gr.add_edge(ss, tt);
		gr.set_edge_weight(p, w);
		edge_info ei;
		ei.type = type;
		gr.set_edge_info(p, ei);
	}
	return 0;
}

int sgraph_compare::search_splice_graph(splice_graph &gr, int32_t p)
{
	if(gr.num_vertices() <= 2) return -1;
	int l = 1;
	int r = gr.num_vertices() - 2;
	while(l <= r)
	{
		int m = (l + r) / 2;
		int32_t p1 = (int32_t)(gr.get_vertex_info(m).lpos);
		int32_t p2 = (int32_t)(gr.get_vertex_info(m).rpos);
		assert(p1 < p2);
		if(p >= p1 && p < p2) return m;
		if(p < p1) r = m - 1;
		if(p >= p2) l = m + 1;
	}
	return -1;
}

int sgraph_compare::identify_unique_boundary_edges()
{
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr3.out_edges(0); it1 != it2; it1++)
	{
		bool b = verify_unique_5end_edge(gr3, *it1);
		if(b == false) continue;

		double w = gr3.get_edge_weight(*it1);
		edge_info ei = gr3.get_edge_info(*it1);

		string s = (ei.type == 3) ? "reference" : "prediction";
		printf("unique 5end %s edge, weight = %.2lf\n", s.c_str(), w);
	}

	for(tie(it1, it2) = gr3.in_edges(gr3.num_vertices() - 1); it1 != it2; it1++)
	{
		bool b = verify_unique_3end_edge(gr3, *it1);
		if(b == false) continue;

		double w = gr3.get_edge_weight(*it1);
		edge_info ei = gr3.get_edge_info(*it1);

		string s = (ei.type == 3) ? "reference" : "prediction";
		printf("unique 3end %s edge, weight = %.2lf\n", s.c_str(), w);
	}

	return 0;
}

bool sgraph_compare::verify_unique_5end_edge(splice_graph &gr, edge_descriptor e)
{
	int t = e->target();
	int type = gr.get_edge_info(e).type;

	// left extend
	int32_t p = gr.get_vertex_info(t).lpos;
	for(int i = t; i >= 1; i--)
	{
		vertex_info vi = gr.get_vertex_info(i);
		int32_t pp = gr.get_vertex_info(i).rpos;
		if(i < t && pp != p) break;

		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.in_edges(i); it1 != it2; it1++)
		{
			if((*it1)->source() != 0) continue;
			int type2 = gr.get_edge_info(*it1).type;
			if(type2 != type) return false;
		}

		p = gr.get_vertex_info(i).lpos;
	}

	// right extend
	p = gr.get_vertex_info(t).rpos;
	for(int i = t + 1; i < gr.num_vertices() - 1; i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		int32_t pp = gr.get_vertex_info(i).lpos;
		if(pp != p) break;

		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.in_edges(i); it1 != it2; it1++)
		{
			if((*it1)->source() != 0) continue;
			int type2 = gr.get_edge_info(*it1).type;
			if(type2 != type) return false;
		}

		p = gr.get_vertex_info(i).rpos;
	}

	return true;
}

bool sgraph_compare::verify_unique_3end_edge(splice_graph &gr, edge_descriptor e)
{
	int s = e->source();
	int type = gr.get_edge_info(e).type;

	// left extend
	int32_t p = gr.get_vertex_info(s).lpos;
	for(int i = s; i >= 1; i--)
	{
		vertex_info vi = gr.get_vertex_info(i);
		int32_t pp = gr.get_vertex_info(i).rpos;
		if(i < s && pp != p) break;

		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.out_edges(i); it1 != it2; it1++)
		{
			if((*it1)->target() != gr.num_vertices() - 1) continue;
			int type2 = gr.get_edge_info(*it1).type;
			if(type2 != type) return false;
		}

		p = gr.get_vertex_info(i).lpos;
	}

	// right extend
	p = gr.get_vertex_info(s).rpos;
	for(int i = s + 1; i < gr.num_vertices() - 1; i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		int32_t pp = gr.get_vertex_info(i).lpos;
		if(pp != p) break;

		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.out_edges(i); it1 != it2; it1++)
		{
			if((*it1)->target() != gr.num_vertices() - 1) continue;
			int type2 = gr.get_edge_info(*it1).type;
			if(type2 != type) return false;
		}

		p = gr.get_vertex_info(i).rpos;
	}

	return true;
}

int sgraph_compare::draw(splice_graph &gr, const string &file)
{
	ofstream fout(file.c_str());
	if(fout.fail())
	{
		printf("open file %s error.\n", file.c_str());
		return 0;
	}

	draw_header(fout);

	double len = 4.0;

	fout<<"\\def\\len{"<<len<<"cm}\n";

	// draw file name
	fout<<"\\node[draw, thick, red] at (1.6 * \\len, 0.58 * \\len) {"<<file.c_str()<<"};\n";

	// draw vertices
	char sx[1024];
	char sy[1024];
	double pos = 0;
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		int d = gr.degree(i);
		//if(d == 0) continue;

		vertex_info vi = gr.get_vertex_info(i);

		pos++;

		sprintf(sx, "s%d", i);
		string s = "";
		fout.precision(0);
		fout<<fixed;
		fout<<"\\node[mycircle, \\colx, draw, label = below:{";
		fout<<vi.lpos % 100000<<"-"<<vi.rpos % 100000;
		fout<<"}] ("<<sx<<") at ("<<pos<<" *\\len, 0.0) {"<<i<<"};\n";
	}

	// draw reference edges
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		edge_info ei = gr.get_edge_info(e);

		int d = (int)(gr.get_edge_weight(e) * 2.0);

		if(ei.type == 2 || ei.type == 4) continue;

		int s = e->source();
		int t = e->target();

		sprintf(sx, "s%d", s);
		sprintf(sy, "s%d", t);

		double bend = 0;
		if(ei.type == 3) bend = -40;

		string line = "line width = 0.12cm, gray, ";

		fout<<"\\draw[->,"<< line.c_str() <<"bend right = "<< bend <<"] ("<<sx<<") to node[gray, label=below:{" << d << "}]{} "<<"("<<sy<<");\n";
	}

	// draw evaluated edges
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		edge_info ei = gr.get_edge_info(e);

		int d = (int)(gr.get_edge_weight(e));

		if(ei.type == 1 || ei.type == 3) continue;

		int s = e->source();
		int t = e->target();

		sprintf(sx, "s%d", s);
		sprintf(sy, "s%d", t);

		double bend = 0;
		if(ei.type == 4) bend = -40;

		string line = "line width = 0.02cm, red,";

		fout<<"\\draw[->,"<< line.c_str() <<"bend right = "<< bend <<"] ("<<sx<<") to node[red, label=above:{" << d << "}]{} "<<"("<<sy<<");\n";
	}

	draw_footer(fout);

	fout.close();
	return 0;
}
