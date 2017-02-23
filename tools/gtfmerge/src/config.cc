#include "config.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

int algo = 2;
double min_transcript_length = 0;
bool multiple_exon = true;

int parse_parameters(int argc, const char ** argv)
{
	for(int i = 1; i < argc; i++)
	{
		if(string(argv[i]) == "-a")
		{
			algo = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-l")
		{
			min_transcript_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-s")
		{
			multiple_exon = false;
		}
	}
	return 0;
}
