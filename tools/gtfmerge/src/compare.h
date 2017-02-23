#ifndef __COMPARE_H__
#define __COMPARE_H__

#include "item.h"
#include "transcript.h"
#include "gene.h"
#include "genome.h"

using namespace std;

typedef pair<transcript, bool> PTB;
typedef pair<double, double> PDD;

// compare two transcripts
bool compare_structure(const transcript &x, const transcript &y);
bool compare_expression(const transcript &x, const transcript &y);
bool compare_intron_chain(const transcript &x, const transcript &y);

// comapre two genes
int compare_transcripts(const vector<transcript> &x, const vector<transcript> &y);
int compare_transcripts(const vector<transcript> &x, const vector<transcript> &y, vector<bool> &v);

// compare two genomes~(two gtf files)
int compare_genome1(const genome &x, const genome &y);
int compare_genome2(const genome &x, const genome &y);

// prepare data for ROC
int ROC(vector<PTB> &vv, int xtotal);
bool transcript_cmp(const PTB &x, const PTB &y);

#endif
