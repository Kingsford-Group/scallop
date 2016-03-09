#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include "manager.h"
#include "config.h"
#include "splice_graph.h"

using namespace std;

int main(int argc, char **argv)
{
	if(argc < 3) return 0;

	load_config(argv[1]);

	manager sc;
	sc.process(argv[2]);

    return 0;
}
