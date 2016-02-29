#ifndef __DECOMPOSER_H__
#define __DECOMPOSER_H__

#include <vector>

using namespace std;

class decomposer 
{
public:
	decomposer(const vector<int> &s, const vector<int> &t);

public:
	// given two sets
	const vector<int> &s;
	const vector<int> &t;
};

#endif
