#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "scallop.h"
#include "config.h"

#include "fheap.h"

using namespace std;

int main(int argc, char **argv)
{
	test_fheap();
	return 0;

	if(argc != 3) return 0;

	load_config(argv[1]);

	scallop sc;
	sc.process(argv[2]);

    return 0;
}
