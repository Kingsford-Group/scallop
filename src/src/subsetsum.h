/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __SUBSETSUM4_H__
#define __SUBSETSUM4_H__

#include <vector>
#include "equation.h"

using namespace std;

typedef pair<int, int> PI;

// partition s and t into s1/s2 and t1/t2
// such that sum(s1) is close to sum(t1)
// AND sum(s2) is close to sum(t2)
class subsetsum
{
public:
	subsetsum(const vector<PI> &s, const vector<PI> &t);

private:
	vector<PI> source;					// given input
	vector<PI> target;					// given target numbers
	int ubound1;						// ubound for source
	int ubound2;						// ubound for target
	vector< vector<int> > table1;		// dp table2
	vector< vector<int> > table2;		// dp table2

public:
	equation eqn;

public:
	int solve();
	int print();
	static int test();

private:
	int rescale();
	int init(const vector<PI> &vv, vector< vector<int> > &table, int ubound);
	int fill(const vector<PI> &vv, vector< vector<int> > &table, int ubound);
	int backtrace(int vi, const vector<PI> &vv, const vector< vector<int> > &table, vector<int> &ss);
	int optimize();
};

#endif
