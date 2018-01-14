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
#include "param_advising.h"

using namespace std;

int main(int argc, const char **argv)
{
	srand(time(0));

	if(argc == 1)
	{
		(new config())->print_copyright();
		print_help();
		printf("\n");
		print_logo();
		return 0;
	}
	vector<config> cfg(1);
	cfg[0].parse_arguments(argc, argv);

	if(cfg[0].verbose >= 1)
	{
		cfg[0].print_copyright();
		printf("\n");
		cfg[0].print_command_line(argc, argv);
		printf("\n");
		//print_parameters();
	}

	if(config::library_type == EMPTY || config::preview_only == true)
	{
		previewer pv(cfg[0]);
		pv.preview();
	}

	//cfg.push_back(config(cfg[0]));
  //cfg[1].update_from_file(string("testing_data/SRR387661_tophat.config").c_str());
  //cfg[0].update_from_file(string("default.config").c_str());
	//cfg.push_back(config(cfg[0]));
	//cfg.push_back(config(cfg[0]));
	//cfg.push_back(config(cfg[0]));

	if(config::preview_only == true) return 0;

	assembler* asmb = new assembler(cfg[0]);
	//(asmb->solve(cfg[0]))->write("testing_data/default.gtf");
	//(asmb->solve(cfg[1]))->write("testing_data/best.gtf");
	//assembler best = asmb.solve(cfg[0]);
	asmb->assemble();
	//assembler best = parameter_advising<assembler,assembler,config>(asmb,cfg);

  //config::output_file = "none.gff";
	//best.write("testing_data/default.gtf");
  //best.write(config::output_file.c_str());
  asmb->write(config::output_file.c_str());

	return 0;
}
