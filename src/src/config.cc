#include "config.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

//// parameters

// for bam file
int32_t min_bundle_gap = 50;
int min_num_hits_in_bundle = 20;
int32_t min_splice_boundary_hits = 1;
uint32_t min_max_splice_boundary_qual = 3;
double min_average_overlap = 2;
int min_max_region_overlap = 5;
double min_region_coverage = 0.5;
int max_num_bundles = -1;
int slope_bin_size = 10;
int slope_min_bin_num = 9;
int slope_min_distance = 100;
int min_slope_score = 45;
int32_t average_read_length = 100;
int slope_std_bin_num = 30;
int32_t average_slope_length = 300;
int pseudo_length_count = 10;
double min_boundary_edge_weight_ratio = 0.05;

// for algorithm
int max_dp_table_size = 10000;
int max_num_subsetsum_solutions = 10;
double max_equation_error_ratio = 0.3;

//// from command line
string algo;
string input_file;
string output_file;

bool output_tex_files;
string fixed_gene_name;
int min_gtf_transcripts_num;
bool fast_mode;

// for simulation
int simulation_num_vertices;
int simulation_num_edges;
int simulation_max_edge_weight;

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
			output_file = string(argv[i + 1]);
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
		else if(string(argv[i]) == "-f")
		{
			fast_mode = true;
		}
		else if(string(argv[i]) == "-s")
		{
			min_gtf_transcripts_num = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-sv")
		{
			simulation_num_vertices = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-se")
		{
			simulation_num_edges = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-sw")
		{
			simulation_max_edge_weight = atoi(argv[i + 1]);
			i++;
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
		else if(strcmp(key, "min_max_splice_boundary_qual")==0)
		{
			min_max_splice_boundary_qual = (uint32_t)atoi(value);
		}
		else if(strcmp(key, "average_read_length")==0)
		{
			average_read_length = (int32_t)atoi(value);
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
		else if(strcmp(key, "max_num_subsetsum_solutions")==0)
		{
			max_num_subsetsum_solutions = (int)atoi(value);
		}
		else if(strcmp(key, "max_equation_error_ratio")==0)
		{
			max_equation_error_ratio = atof(value);
		}
		else if(strcmp(key, "slope_bin_size")==0)
		{
			slope_bin_size = (int)atoi(value);
		}
		else if(strcmp(key, "slope_min_bin_num")==0)
		{
			slope_min_bin_num = (int)atoi(value);
		}
		else if(strcmp(key, "slope_std_bin_num")==0)
		{
			slope_std_bin_num = (int)atoi(value);
		}
		else if(strcmp(key, "slope_min_distance")==0)
		{
			slope_min_distance = (int)atoi(value);
		}
	}

	return 0;
}
