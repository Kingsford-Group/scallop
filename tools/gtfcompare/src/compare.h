#ifndef __COMPARE_H__
#define __COMPARE_H__

#include "exon.h"
#include "transcript.h"
#include "gene.h"
#include "genome.h"

using namespace std;

// compare two transcripts
bool compare_transcript_structure(const transcript &x, const transcript &y);
bool compare_transcript_expression(const transcript &x, const transcript &y);
bool compare_transcript(const transcript &x, const transcript &y);

// comapre two genes
int compare_gene(const gene &x, const gene &y);

// compare two genomes~(two gtf files)
int compare_genome(const genome &x, const genome &y);

#endif
