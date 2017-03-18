/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __GRAPH_BASE_H__
#define __GRAPH_BASE_H__

#include <vector>
#include <map>
#include <string>

#include "vertex_base.h"
#include "edge_base.h"

using namespace std;

typedef map<int, string> MIS;
typedef pair<int, string> PIS;
typedef map<edge_descriptor, string> MES;
typedef pair<edge_descriptor, string> PES;
typedef map<edge_descriptor, bool> MEB;
typedef map<edge_descriptor, double> MED;
typedef pair<edge_descriptor, double> PED;
typedef map<edge_descriptor, int> MEI;
typedef pair<edge_descriptor, int> PEI;
typedef vector<edge_descriptor> VE;
typedef set<edge_descriptor> SE;

class graph_base
{
public:
	graph_base();
	graph_base(const graph_base &gr);
	virtual ~graph_base();

protected:
	vector<vertex_base*> vv;
	set<edge_base*> se;

public:
	// modify the graph
	virtual int copy(const graph_base &gr);
	virtual int add_vertex();
	virtual int clear_vertex(int v);
	virtual int clear();
	virtual edge_descriptor add_edge(int s, int t) = 0;
	virtual int remove_edge(edge_descriptor e) = 0;
	virtual int remove_edge(int s, int t) = 0;

	// access functions
	virtual size_t support_size() const;
	virtual size_t num_vertices() const;
	virtual size_t num_edges() const;
	virtual int degree(int v) const;
	virtual PEB edge(int s, int t);
	virtual PEEI edges() const;
	virtual vector<edge_descriptor> edges(int x, int y);
	virtual set<int> adjacent_vertices(int v);
	virtual PEEI out_edges(int x) = 0;
	virtual int get_edge_indices(VE &e2i, MEI &i2e);

	// algorithms
	virtual int bfs(int s, vector<int> &v);
	virtual int bfs(int s, vector<int> &v, vector<int> &b);
	virtual int bfs(int s, set<edge_descriptor> &ss);
	virtual bool bfs(const vector<int> &vs, int t, const set<edge_descriptor> &fb);
	virtual bool check_path(int s, int t);
	virtual bool compute_shortest_path(int s, int t, vector<int> &p);
	virtual bool check_nested();
	virtual bool intersect(edge_descriptor ex, edge_descriptor ey) = 0;

	// draw
	virtual int draw(const string &file, const MIS &mis, const MES &mes, double len) = 0;
	virtual int print() const;
};

#endif
