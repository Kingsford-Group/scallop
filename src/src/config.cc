#include "config.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

int min_num_hits_in_bundle;
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
int min_max_region_overlap;
double min_region_coverage;
int max_num_bundles;
int max_dp_table_size;

string algo;
string input_file;
string output_gtf_file;
bool output_tex_files;
string fixed_gene_name;

bool parse_arguments(int argc, const char ** argv)
{
	output_tex_files = false;
	bool b = false;
	for(int i = 1; i < argc; i++)
	{
		if(string(argv[i]) == "-c")
		{
			load_config(argv[i + 1]);
			b = true;
			i++;
		}
		else if(string(argv[i]) == "-a")
		{
			algo = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-o")
		{
			output_gtf_file = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-i")
		{
			input_file = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-g")
		{
			fixed_gene_name = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-t")
		{
			output_tex_files = true;
		}
	}

	return b;
}

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

		if(strcmp(key, "min_num_hits_in_bundle")==0)
		{
			min_num_hits_in_bundle = (int)atoi(value);
		}
		else if(strcmp(key, "min_bundle_gap")==0)
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
			min_average_overlap = (double)atof(value);
		}
		else if(strcmp(key, "min_max_region_overlap")==0)
		{
			min_max_region_overlap = (int)atoi(value);
		}
		else if(strcmp(key, "min_region_coverage")==0)
		{
			min_region_coverage = (double)atof(value);
		}
		else if(strcmp(key, "max_num_bundles")==0)
		{
			max_num_bundles = (int)atoi(value);
		}
		else if(strcmp(key, "max_dp_table_size")==0)
		{
			max_dp_table_size = (int)atoi(value);
		}
	}

	return 0;
}
