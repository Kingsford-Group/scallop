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
int32_t min_left_boundary_hits;
int32_t min_right_boundary_hits;
uint32_t min_max_splice_boundary_qual;
uint32_t min_max_left_boundary_qual;
uint32_t min_max_right_boundary_qual;
int32_t average_read_length;
uint32_t min_boundary_score;
int32_t ascending_step;
int32_t descending_step;
uint32_t min_ascending_score;
uint32_t min_descending_score;
int num_sample_positions;
double min_average_overlap;

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
		else if(strcmp(key, "min_left_boundary_hits")==0)
		{
			min_left_boundary_hits = (int32_t)atoi(value);
		}
		else if(strcmp(key, "min_right_boundary_hits")==0)
		{
			min_right_boundary_hits = (int32_t)atoi(value);
		}
		else if(strcmp(key, "min_max_splice_boundary_qual")==0)
		{
			min_max_splice_boundary_qual = (uint32_t)atoi(value);
		}
		else if(strcmp(key, "min_max_left_boundary_qual")==0)
		{
			min_max_left_boundary_qual = (uint32_t)atoi(value);
		}
		else if(strcmp(key, "min_max_right_boundary_qual")==0)
		{
			min_max_right_boundary_qual = (uint32_t)atoi(value);
		}
		else if(strcmp(key, "average_read_length")==0)
		{
			average_read_length = (int32_t)atoi(value);
		}
		else if(strcmp(key, "ascending_step")==0)
		{
			ascending_step = (int32_t)atoi(value);
		}
		else if(strcmp(key, "descending_step")==0)
		{
			descending_step = (int32_t)atoi(value);
		}
		else if(strcmp(key, "min_boundary_score")==0)
		{
			min_boundary_score = (uint32_t)atoi(value);
		}
		else if(strcmp(key, "min_ascending_score")==0)
		{
			min_ascending_score = (uint32_t)atoi(value);
		}
		else if(strcmp(key, "min_descending_score")==0)
		{
			min_descending_score = (uint32_t)atoi(value);
		}
		else if(strcmp(key, "num_sample_positions")==0)
		{
			num_sample_positions = (int)atoi(value);
		}
		else if(strcmp(key, "min_average_overlap")==0)
		{
			min_average_overlap = (int)atof(value);
		}

	}

	return 0;
}
