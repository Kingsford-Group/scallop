/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "util.h"
#include <stdint.h>
#include <map>
#include <sstream>

using namespace std;

//// constants
#define MAX_NUM_CIGAR 7

#define START_BOUNDARY 1
#define END_BOUNDARY 2
#define LEFT_SPLICE 3
#define RIGHT_SPLICE 4
#define LEFT_RIGHT_SPLICE 5
#define MIDDLE_CUT 6

#define TRIVIAL 0
#define NORMAL 1

// five types for decomposition
#define SMALLEST_EDGE 0
#define NEGLIGIBLE_EDGE 1
#define SPLITTABLE_SIMPLE 2
#define SPLITTABLE_HYPER 3
#define UNSPLITTABLE_SINGLE 4
#define UNSPLITTABLE_MULTIPLE 5
#define TRIVIAL_VERTEX 6

#define EMPTY -1
#define UNSTRANDED 0
#define FR_FIRST 1
#define FR_SECOND 2

#define BRIDGE -18

class config{
public:
  config();

  //// parameters
  // for bam file and reads
  int min_flank_length;
  int max_edit_distance;
  int32_t min_bundle_gap;
  int min_num_hits_in_bundle;
  uint32_t min_mapping_quality;
  int32_t min_splice_boundary_hits;
  bool uniquely_mapped_only;
  bool use_second_alignment;

  // for preview
  bool preview_only;
  int max_preview_reads;
  int max_preview_spliced_reads;
  int min_preview_spliced_reads;
  double preview_infer_ratio;

  // for identifying subgraphs
  int32_t min_subregion_gap;
  double min_subregion_overlap;
  int32_t min_subregion_length;
  int min_subregion_ladders;

  // for subsetsum and router
  int max_dp_table_size;
  int min_router_count;

  // for splice graph
  double max_intron_contamination_coverage;
  double min_surviving_edge_weight;
  double max_decompose_error_ratio[7];
  double min_transcript_numreads;
  double min_transcript_coverage;
  double min_single_exon_coverage;
  double min_transcript_coverage_ratio;
  int min_transcript_length_base;
  int min_transcript_length_increase;
  int min_exon_length;
  int max_num_exons;

  // for simulation
  int simulation_num_vertices;
  int simulation_num_edges;
  int simulation_max_edge_weight;

  // input and output
  string algo;
  string input_file;
  string ref_file;
  string ref_file1;
  string ref_file2;
  string output_file;

  // for controling
  bool output_tex_files;
  string fixed_gene_name;
  int max_num_bundles;
  int library_type;
  int min_gtf_transcripts_num;
  int batch_bundle_size;
  int verbose;
  string version;

  // parse arguments
  int print_command_line(int argc, const char ** argv);
  int parse_arguments(int argc, const char ** argv);
  int print_parameters();
  int print_copyright();
  void update_from_file(char * fname);
};

int print_logo();
int print_help();

#endif
