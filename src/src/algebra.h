#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__

#include "splice_graph.h"

#define SMIN 0.0000001

typedef vector< vector<double> > VVI;

class algebra
{
public:
	algebra();
	~algebra();

public:
	splice_graph gr;	// splice graph

public:
	VE i2e;
	MEI e2i;
	VVI pem0;			// backup
	VVI pem;			// path-edge-matrxi
	int rank;			// rank of pem

public:
	int solve(int n, int m);
	int solve(const string &file);
	int process();
	int simulate(int n, int m);
	int build_path_edge_matrix();
	int build_path_edge_matrix(int s, vector<int> &v);
	int sample_valid_paths();
	int sample_paths();
	int gaussian_elimination();
	int gaussian_elimination(int r);
	int choose_main_element(int r);
	int choose_main_row(int r);
	int choose_main_column(int r);
	int normalize_row(int r);
	int eliminate_column(int r);
	int add_to_row(int s, int t, double c);
	int exchange_row(int s, int t);
	int exchange_column(int s, int t);
	int check_edges_cover();
	int check();

public:
	int print();
};

#endif
