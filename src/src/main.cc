#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include "assembler.h"
#include "config.h"
#include "subsetsum.h"

using namespace std;

int main(int argc, const char **argv)
{
	bool b = parse_arguments(argc, argv);
	if(b == false) return 0;

	assembler asmbl;
	asmbl.process();

    return 0;
}
