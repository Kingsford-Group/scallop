#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "manager.h"
#include "config.h"
#include "subsetsum.h"
#include "path_space.h"

using namespace std;

int main(int argc, char **argv)
{
	if(argc < 3) return 0;
	path_space ps;
	ps.solve(atoi(argv[1]), atoi(argv[2]));
	return 0;

	load_config(argv[1]);

	manager sc;
	sc.process(argv[2]);

    return 0;
}
