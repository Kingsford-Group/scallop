#ifndef __ASSEMBER_BASE_H__
#define __ASSEMBER_BASE_H__

#include "splice_graph.h"
#include "path.h"

class assember
{
public:
	assember(splice_graph &g);
	virtual ~assember();

public:
	splice_graph &gr;				// splice graph
	vector<path> paths0;			// paths from greedy algorithm
	vector<path> paths1;			// paths from (our) iterate algorithm

public:
	int solve();
	int draw(const string &file) const;
	int print() const;

private:
	int greedy();
	int iterate();

	int update_weights();

	path compute_maximum_forward_path() const;
	path compute_maximum_path() const;
	int decrease_path(const path &p);
	int increase_path(const path &p);
	int add_backward_path(const path &p);
	int remove_backward_path(const path &p);
	int resolve(const path &px, const path &py, path &qx, path &qy) const;
	double compute_bottleneck_weight(const path &p) const;

	int backup_edge_weights(MED &med) const;
	int recover_edge_weights(const MED &med);
};

#endif
