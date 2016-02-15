#ifndef __SGRAPH_H__
#define __SGRAPH_H__

#include "bundle.h"
#include "dgraph.h"
#include "path.h"

// class for splice graph
class sgraph
{
public:
	sgraph();
	virtual ~sgraph();

public:
	dgraph gr;						// splice graph
	vector<path> paths0;			// paths from greedy algorithm
	vector<path> paths1;			// paths from (our) iterate algorithm

public:
	int solve();
	int load(const string &file);
	int build(const bundle &bd);
	int draw(const string &file);
	int print();

private:
	int update_weights();

	int greedy();
	int iterate();

	path compute_maximum_forward_path();
	path compute_maximum_path();
	int decrease_path(const path &p);
	int increase_path(const path &p);
	int add_backward_path(const path &p);
	int remove_backward_path(const path &p);
	int resolve(const path &px, const path &py, path &qx, path &qy) const;
	double compute_bottleneck_weight(const path &p);

	int backup_edge_weights(MED &med);
	int recover_edge_weights(const MED &med);
};

#endif
