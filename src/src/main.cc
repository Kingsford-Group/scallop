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
#include "interval_map.h"

using namespace std;

int main(int argc, const char **argv)
{
	srand(time(0));

	if(argc == 1)
	{
		print_help();
		printf("\n");
		print_logo();
		return 0;
	}

	print_command_line(argc, argv);
	parse_arguments(argc, argv);
	print_parameters();

	if(input_file == "") return 0;

	string s = input_file.substr(input_file.size() - 3, 3);
	if(s != "bam") return 0;

	assembler asmb;
	asmb.assemble();

	return 0;
}

