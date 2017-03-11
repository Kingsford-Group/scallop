#ifndef __SCALLOP3_H__
#define __SCALLOP3_H__

#include "splice_graph.h"
#include "hyper_set.h"
#include "equation.h"
#include "router.h"
#include "path.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair< vector<int>, vector<int> > PVV;
typedef pair<PEE, int> PPEEI;
typedef map<PEE, int> MPEEI;
typedef pair<int, int> PI;
typedef map<int, int> MI;

// for noisy splice graph
class scallop
{
public:
	scallop();
	scallop(const string &name, const splice_graph &gr, const hyper_set &hs);
	virtual ~scallop();

public:
	int assemble();

public:
	string name;						// name for this gene
	splice_graph gr;					// splice graph
	MEI e2i;							// edge map, from edge to index
	VE i2e;								// edge map, from index to edge
	MEV mev;							// super edges
	vector<int> v2v;					// vertex map
	hyper_set hs;						// hyper edges
	int round;							// iteration
	set<int> nonzeroset;				// vertices with degree >= 1
	vector<path> paths;					// predicted transcripts

private:
	// init
	int classify();
	int init_vertex_map();
	int init_super_edges();
	int init_inner_weights();
	int init_nonzeroset();
	int add_pseudo_hyper_edges();
	int refine_splice_graph();

	// resolve iteratively
	bool resolve_trivial_vertex(int type, double jump_ratio);
	bool resolve_trivial_vertex_fast(double jump_ratio);
	bool resolve_single_trivial_vertex_fast(int i, double jump_ratio);
	bool resolve_small_edges(double max_ratio);
	bool resolve_splittable_vertex(int type, int degree, double max_ratio);
	bool resolve_unsplittable_vertex(int type, int degree, double max_ratio);
	bool resolve_hyper_edge(int fsize);

	// smooth vertex
	int balance_vertex(int x);
	double compute_balance_ratio(int x);

	// decomposing subroutines
	int compute_smallest_edge(int x, double &ratio);
	int decompose_trivial_vertex(int v);
	int decompose_vertex_extend(int v, MPID &pe2w);
	int decompose_vertex_replace(int v, MPID &pe2w);
	int classify_trivial_vertex(int v, bool fast);
	int exchange_sink(int old_sink, int new_sink);
	int split_vertex(int x, const vector<int> &xe, const vector<int> &ye);
	int split_edge(int exi, double w);
	int merge_adjacent_edges(int x, int y);
	int merge_adjacent_edges(int x, int y, double ww);
	int merge_adjacent_equal_edges(int x, int y);
	int remove_edge(int e);
	int split_merge_path(const VE &p, double w);
	int split_merge_path(const vector<int> &p, double w);
	int collect_existing_st_paths();
	int collect_path(int e);
	int compute_length(const path &p);
	int greedy_decompose();

	// stats, print, and draw
	int print();
	int stats();
	int summarize_vertices();
	int draw_splice_graph(const string &file);
	vector<int> topological_sort();
};

#endif
