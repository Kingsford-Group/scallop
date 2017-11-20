/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __HYPER_SET_H__
#define __HYPER_SET_H__

#include <map>
#include <set>
#include <vector>

#include "util.h"
#include "directed_graph.h"

using namespace std;

typedef pair<vector<int>, int> PVII;
typedef map<vector<int>, int> MVII;
typedef map< int, set<int> > MISI;
typedef pair< int, set<int> > PISI;
typedef vector< vector<int> > VVI;
typedef map< pair<int, int>, int> MPII;
typedef pair< pair<int, int>, int> PPII;

class hyper_set
{
public:
	MVII nodes;			// hyper-edges using list-of-nodes
	VVI edges;			// hyper-edges using list-of-edges
	vector<int> ecnts;	// counts for edges
	MISI e2s;			// index: from edge to hyper-edges

public:
	int clear();
	int add_node_list(const set<int> &s);
	int add_node_list(const set<int> &s, int c);
	int add_node_list(const vector<int> &s, int c);
	int build(directed_graph &gr, MEI &e2i);
	int build_edges(directed_graph &gr, MEI &e2i);
	int build_index();
	int update_index();
	set<int> get_intersection(const vector<int> &v);
	MI get_successors(int e);
	MI get_predecessors(int e);
	MPII get_routes(int x, directed_graph &gr, MEI &e2i);
	int print();

public:
	int replace(int x, int e);
	int replace(int x, int y, int e);
	int replace(const vector<int> &x, int e);
	int remove(int e);
	int remove(const vector<int> &x);
	int remove(const set<int> &x);
	int remove_pair(int x, int y);
	int insert_between(int x, int y, int e);
	bool useful(const vector<int> &v, int k1, int k2);
	bool extend(int e);
	bool left_extend(int e);
	bool left_extend(const vector<int> &s);
	bool right_extend(int e);
	bool right_extend(const vector<int> &s);
	bool left_dominate(int e);
	bool right_dominate(int e);
};

#endif
