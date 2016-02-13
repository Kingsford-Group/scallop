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
	int print(int index);
	int clear();
};

#endif
