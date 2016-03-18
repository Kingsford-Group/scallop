#ifndef __GRAPH_BASE_H__
#define __GRAPH_BASE_H__

#include <vector>
#include <map>

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
	virtual int add_vertex();
	virtual edge_descriptor add_edge(int s, int t);
	virtual int remove_edge(edge_descriptor e);
	virtual int clear_vertex(int v);
	virtual int clear();
	virtual int remove_edge(int s, int t) = 0;

	// access functions
	virtual size_t num_vertices() const;
	virtual size_t num_edges() const;
	virtual int degree(int v) const;
	virtual PEE edges() const;
	virtual PEE out_edges(int x) const = 0;
	virtual PEB edge(int s, int t) const = 0;
	virtual vector<edge_descriptor> edges(int x, int y) const = 0;
	virtual set<int> adjacent_vertices(int v) const = 0;

	// algorithms
	virtual int bfs(int s, vector<int> &v) const;
	virtual int bfs(int s, vector<int> &v, vector<int> &b) const;
	virtual bool check_path(int s, int t) const;
	virtual bool compute_shortest_path(int s, int t, vector<int> &p) const;
	virtual bool check_nested() const;
	virtual bool intersect(edge_descriptor ex, edge_descriptor ey) const = 0;

	// draw
	virtual int draw(const string &file, const MIS &mis, const MES &mes, double len) const = 0;
};

#endif
