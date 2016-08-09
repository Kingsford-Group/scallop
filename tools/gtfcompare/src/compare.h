#ifndef __COMPARE_H__
#define __COMPARE_H__

#include "exon.h"
#include "transcript.h"
#include "gene.h"
#include "genome.h"

using namespace std;

// compare two transcripts
bool compare_structure(const transcript &x, const transcript &y);
bool compare_intron_chain(const transcript &x, const transcript &y);
bool compare_expression(const transcript &x, const transcript &y);

// comapre two genes
int compare_gene(const gene &x, const gene &y, int mode);
int compare_transcripts(const vector<transcript> &x, const vector<transcript> &y, int mode);

// compare two genomes~(two gtf files)
int compare_genome1(const genome &x, const genome &y);
int compare_genome2(const genome &x, const genome &y);

// remove single exon transcripts
int remove_single_exon_transcripts(genome &gm);

#endif
