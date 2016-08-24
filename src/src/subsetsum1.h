#ifndef __SUBSETSUM1_H__
#define __SUBSETSUM1_H__

#include <vector>
#include "equation.h"

using namespace std;

typedef pair<int, int> PI;

// for given s compute all subsetsums 
// and then backtrace for each given t
class subsetsum1
{
public:
	subsetsum1(const vector<PI> &s, const vector<PI> &t);

private:
	vector<PI> source;					// given input
	vector<PI> target;					// given target numbers
	int ubound;							// upper bound we would try 

	vector< vector<bool> > table1;		// dp table1
	vector< vector<int> > table2;		// dp table2

public:
	vector<equation> eqns;

public:
	int solve();
	int print();
	static int test();

private:
	int init();
	int fill();
	int optimize();
	int backtrace(int ti, equation &eqn);
};

#endif
