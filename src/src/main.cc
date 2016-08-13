#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <ctime>

#include "assembler.h"
#include "config.h"

using namespace std;

int main(int argc, const char **argv)
{
	srand(time(0));
	parse_arguments(argc, argv);
	print_parameters();


	assembler asmbl;
	asmbl.process();

    return 0;
}
