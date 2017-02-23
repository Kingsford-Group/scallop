#ifndef __CUFFROC_H__
#define __CUFFROC_H__

#include "cuffitem.h"
#include "quantitem.h"
#include "genome.h"
#include <vector>
#include <map>

using namespace std;

class gtfcuff
{
public:
	gtfcuff(const string &cufffile);

public:
	vector<cuffitem> items;
	vector<transcript> vpred;	// predicted transcripts
	vector<transcript> vref;	// reference transcripts
	vector<quantitem> qitems;	// quant items (for Salmon)
	map<string, int> t2i;		// transcript_id to item-index
	map<string, int> t2p;		// transcript_id to pred-index
	map<string, int> t2r;		// transcript_id to ref-index
	map<string, int> t2q;		// transcript_id to quant-index

public:
	int read_cuff(const string &file);
	int read_quant(const string &file);
	int assign_ref(const string &file);
	int assign_pred(const string &file);
	int build_cuff_index();
	int build_pred_index();
	int build_ref_index();
	int build_quant_index();

	int roc(int refsize);
	int acc(int refsize);
	int roc_trunc(int refsize, double min_coverage, double max_coverage);
	int roc_quant(const string &qfile, double min_tpm, double max_tpm);
	int acc_quant(const string &qfile, double tpm_threshold);
	int split(const string &f1, const string &f2);
	int classify();
	int quant();
};

#endif
