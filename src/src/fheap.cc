#include "fheap.h"
#include <cstdio>

int test_fheap()
{
	fheap f;
	handle_t h;
	vector<handle_t> v;
	h = f.push(0); v.push_back(h);
	h = f.push(1); v.push_back(h);
	h = f.push(2); v.push_back(h); 
	h = f.push(3); v.push_back(h);
	h = f.push(4); v.push_back(h);

	f.decrease(v[0], 8);
	f.increase(v[4], 2);

	for(int i = 0; i < v.size(); i++)
	{
		printf("%.0lf, ", *(v[i]));
	}
	printf("\n");

	while(f.empty() == false)
	{
		printf("%.0lf, ", f.top());
		f.pop();
	}
	printf("\n");

	return 0;
}
