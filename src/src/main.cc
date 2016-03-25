#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include "manager.h"
#include "config.h"
#include "subsetsum.h"

using namespace std;

int main(int argc, char **argv)
{
	if(argc < 3) return 0;

	load_config(argv[1]);

	manager sc;
	if(argc == 3) sc.process(argv[2], "");
	else sc.process(argv[2], argv[3]);

    return 0;
}
