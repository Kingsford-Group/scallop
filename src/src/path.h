#ifndef __PATH_H__
#define __PATH_H__

#include <vector>

using namespace std;

class path
{
public:
	path();
	~path();

public:
	vector<int> v;
	double abd;

public:
	int clear();
	int print(int index) const;
	vector<int> index(int n) const;
};

#endif
