#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <ctime>

#include "assembler.h"
#include "config.h"
#include "subsetsum.h"

using namespace std;

int main(int argc, const char **argv)
{
	srand(time(0));
	parse_arguments(argc, argv);

	assembler asmbl;
	asmbl.process();

    return 0;
}
