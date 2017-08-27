/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __SGRAPH_MERGE_H__
#define __SGRAPH_MERGE_H__

#include "splice_graph.h"
#include "hyper_set.h"
#include "interval_map.h"

class sgraph_merge
{
public:
	sgraph_merge(const splice_graph &g1, const splice_graph &g2, const hyper_set &h1, const hyper_set &h2);

public:
	splice_graph gr1;		
	splice_graph gr2;		
	splice_graph gr3;		

	hyper_set hs1;
	hyper_set hs2;
	hyper_set hs3;

public:
	split_interval_map imap;

public:
	int compare(const string &texfile = "");

private:
	int build_split_interval_map(splice_graph &gr);
	int add_vertices(splice_graph &gr);
	int add_inner_edges(splice_graph &gt, splice_graph &gr, int type);
	int add_existing_edges(splice_graph &gt, splice_graph &gr, int type);
	int search_splice_graph(splice_graph &gr, int32_t p);

	int compare_splice_positions();
	int compare_boundary_edges();
	bool verify_unique_5end_edge(splice_graph &gr, edge_descriptor e);
	bool verify_unique_3end_edge(splice_graph &gr, edge_descriptor e);
	int draw(splice_graph &gr, const string &file);
};

#endif
