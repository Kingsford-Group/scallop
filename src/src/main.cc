#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "manager.h"
#include "config.h"
#include "subsetsum.h"
#include "decomposer.h"
#include "algebra.h"

using namespace std;

int main(int argc, char **argv)
{
	algebra ab;

	if(argc == 2) ab.solve(argv[1]);
	if(argc == 3) ab.solve(atoi(argv[1]), atoi(argv[2]));

	return 0;

	if(argc < 3) return 0;

	load_config(argv[1]);

	manager sc;
	sc.process(argv[2]);

    return 0;
}
