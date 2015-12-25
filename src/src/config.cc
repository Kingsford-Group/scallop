#include "config.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

int32_t min_bundle_gap;
int32_t min_splice_boundary_hits;
int32_t min_start_boundary_hits;
int32_t min_end_boundary_hits;
uint32_t min_max_splice_boundary_qual;
uint32_t min_max_start_boundary_qual;
uint32_t min_max_end_boundary_qual;
int32_t hits_window_size;
double min_boundary_score;

int load_config(const char * conf_file)
{
	ifstream fin(conf_file);
	if(fin.fail())
	{
		cout<<"open config file error "<<conf_file<<endl;
		return -1;
	}

	char buf[1024];
	char key[1024];
	char value[1024];
	
	while(fin.getline(buf, 1024, '\n'))
	{
		stringstream sstr(buf);
		sstr>>key>>value;

		if(strcmp(key, "min_bundle_gap")==0)
		{
			min_bundle_gap = (int32_t)atoi(value);
		}
		else if(strcmp(key, "min_splice_boundary_hits")==0)
		{
			min_splice_boundary_hits = (int32_t)atoi(value);
		}
		else if(strcmp(key, "min_start_boundary_hits")==0)
		{
			min_start_boundary_hits = (int32_t)atoi(value);
		}
		else if(strcmp(key, "min_end_boundary_hits")==0)
		{
			min_end_boundary_hits = (int32_t)atoi(value);
		}
		else if(strcmp(key, "min_max_splice_boundary_qual")==0)
		{
			min_max_splice_boundary_qual = (uint32_t)atoi(value);
		}
		else if(strcmp(key, "min_max_start_boundary_qual")==0)
		{
			min_max_start_boundary_qual = (uint32_t)atoi(value);
		}
		else if(strcmp(key, "min_max_end_boundary_qual")==0)
		{
			min_max_end_boundary_qual = (uint32_t)atoi(value);
		}
		else if(strcmp(key, "hits_window_size")==0)
		{
			hits_window_size = (int32_t)atoi(value);
		}
		else if(strcmp(key, "min_boundary_score")==0)
		{
			min_boundary_score = (int32_t)atoi(value);
		}
	}

	return 0;
}
