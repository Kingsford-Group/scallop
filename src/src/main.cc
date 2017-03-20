/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <ctime>
#include <cassert>
#include <sstream>

#include "config.h"
#include "assembler.h"
#include "matrix.h"

using namespace std;

int main(int argc, const char **argv)
{
	srand(time(0));

	if(argc == 1)
	{
		print_copyright();
		print_help();
		printf("\n");
		print_logo();
		return 0;
	}

	parse_arguments(argc, argv);

	if(input_file == "") return 0;
	if(output_file == "") return 0;

	assembler asmb;
	asmb.assemble();

	return 0;
}
