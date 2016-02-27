#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "manager.h"
#include "config.h"
#include "subsetsum.h"

using namespace std;

int main(int argc, char **argv)
{
	vector<int> s;
	s.push_back(189);
	s.push_back(958);
	s.push_back(1105);

	vector<int> t;
	t.push_back(440);
	t.push_back(637);
	t.push_back(319);
	t.push_back(851);

	subsetsum sss(s, t);
	sss.print();
	return 0;

	if(argc < 3) return 0;

	load_config(argv[1]);

	manager sc;
	sc.process(argv[2]);

    return 0;
}
