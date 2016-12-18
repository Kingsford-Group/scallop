#ifndef __CUFFROC_H__
#define __CUFFROC_H__

#include "cuffitem.h"
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
	map<string, int> t2i;		// transcript_id to item-index
	map<string, int> t2p;		// transcript_id to pred-index
	map<string, int> t2r;		// transcript_id to ref-index

public:
	int read_cuff(const string &file);
	int assign_ref(const string &file);
	int assign_pred(const string &file);
	int build_cuff_index();
	int build_pred_index();
	int build_ref_index();

	int roc(int refsize);
	int classify(const string &f1, const string &f2);
	int quant();
};

#endif
