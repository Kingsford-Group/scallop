#ifndef __SGRAPH_COMPARE_H__
#define __SGRAPH_COMPARE_H__

#include "splice_graph.h"
#include "interval_map.h"

class sgraph_compare
{
public:
	sgraph_compare(const splice_graph &g1, const splice_graph &g2);

public:
	splice_graph gr1;		// reference sgraph
	splice_graph gr2;		// evaluate sgraph
	splice_graph gr3;		// combination

public:
	split_interval_map imap;

public:
	int compare(const string &file);

private:
	int build_split_interval_map(splice_graph &gr);
	int add_vertices(splice_graph &gr);
	int add_inner_edges(splice_graph &gt, splice_graph &gr, int type);
	int add_existing_edges(splice_graph &gt, splice_graph &gr, int type);
	int search_splice_graph(splice_graph &gr, int32_t p);
	int draw(splice_graph &gr, const string &file);
	int identify_unique_boundary_edges();
	bool verify_unique_5end_edge(splice_graph &gr, edge_descriptor e);
	bool verify_unique_3end_edge(splice_graph &gr, edge_descriptor e);
};

#endif
