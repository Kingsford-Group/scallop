#include "fibonacci_heap.h"
#include <cstdio>
#include <vector>

using namespace std;

int test_fibonacci_heap()
{
	fibonacci_heap f;
	handle_t h;
	vector<handle_t> v;
	h = f.push(fnode(0, 0)); v.push_back(h);
	h = f.push(fnode(1, 1)); v.push_back(h);
	h = f.push(fnode(2, 2)); v.push_back(h); 
	h = f.push(fnode(3, 3)); v.push_back(h);
	h = f.push(fnode(4, 4)); v.push_back(h);

	for(int i = 0; i < v.size(); i++)
	{
		printf("%d -> %.0lf\n", (*v[i]).v, (*v[i]).w);
	}
	printf("\n");

	f.increase(v[0], fnode(0, 8));
	f.decrease(v[4], fnode(4, 2));
	f.push(fnode(5, 5));

	for(int i = 0; i < v.size(); i++)
	{
		printf("%d -> %.0lf\n", (*v[i]).v, (*v[i]).w);
	}
	printf("\n");


	while(f.empty() == false)
	{
		fnode fn = f.top();
		printf("%d -> %.0lf\n", fn.v, fn.w);
		f.pop();
	}

	return 0;
}
