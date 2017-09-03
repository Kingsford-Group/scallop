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
#include "previewer.h"
#include "assembler.h"

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

	if(verbose >= 1)
	{
		print_copyright();
		printf("\n");
		print_command_line(argc, argv);
		printf("\n");
		//print_parameters();
	}

	if(library_type == EMPTY || preview_only == true)
	{
		previewer pv;
		pv.preview();
	}

	if(preview_only == true) return 0;

	assembler asmb;
	asmb.assemble();

	return 0;
}
