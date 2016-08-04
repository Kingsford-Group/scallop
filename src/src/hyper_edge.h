#ifndef __HYPER_EDGE_H__
#define __HYPER_EDGE_H__

#include <vector>

using namespace std;

class hyper_edge
{
public:
	hyper_edge(const vector<int> &v);
	hyper_edge(const vector<int> &v, int c);

	bool operator<(const hyper_edge &h) const;

public:
	vector<int> v;		// vertices
	int count;			// # of segments

public:
	int increase();
	int print(int index) const;
};

#endif
