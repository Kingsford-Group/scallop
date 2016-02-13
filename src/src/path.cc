#include "path.h"

#include <cstdio>

path::path()
{
	v.clear();
	abd = 0;
}

int path::clear()
{
	v.clear();
	abd = 0;
	return 0;
}
	
path::~path()
{}

int path::print(int index)
{
	if(v.size() == 0) return 0;
	printf("path %d: abundance = %.2lf, vertices = ", index, abd);
	for(int i = 0; i < v.size() - 1; i++) printf("%d, ", v[i]);
	printf("%d\n", v[v.size() - 1]);
	return 0;
}
