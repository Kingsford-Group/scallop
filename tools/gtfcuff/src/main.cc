#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <cassert>

#include "item.h"
#include "gtfcuff.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc != 4 && argc != 5 && argc != 6)
	{
		cout<<"usage: " <<endl;
		cout<<"       " <<argv[0] << " roc <cuff.tmap> <ref-size>"<<endl;
		cout<<"       " <<argv[0] << " acc <cuff.tmap> <ref-size>"<<endl;
		cout<<"       " <<argv[0] << " roc-trunc <cuff.tmap> <ref-size> <min-coverage> <max-coverage>"<<endl;
		cout<<"       " <<argv[0] << " roc-quant <cuff.tmap> <quant-file> <min-tpm> <max-tpm>"<<endl;
		cout<<"       " <<argv[0] << " acc-quant <cuff.tmap> <quant-file> <tpm-threshold>"<<endl;
		cout<<"       " <<argv[0] << " classify <cuff.tmap> <pred-gtf-file>"<<endl;
		cout<<"       " <<argv[0] << " quant <cuff.tmap> <pred-gtf-file> <ref-gtf-file>"<<endl;
		cout<<"       " <<argv[0] << " split <cuff.tmap> <pred-gtf-file> <true-file> <false-file>"<<endl;
		return 0;
	}

	if(string(argv[1]) == "roc")
	{
		assert(argc == 4);
		gtfcuff cuff(argv[2]);
		cuff.roc(atoi(argv[3]));
	}

	if(string(argv[1]) == "acc")
	{
		assert(argc == 4);
		gtfcuff cuff(argv[2]);
		cuff.acc(atoi(argv[3]));
	}

	if(string(argv[1]) == "roc-trunc")
	{
		assert(argc == 6);
		gtfcuff cuff(argv[2]);
		cuff.roc_trunc(atoi(argv[3]), atof(argv[4]), atof(argv[5]));
	}

	if(string(argv[1]) == "roc-quant")
	{
		assert(argc == 6);
		gtfcuff cuff(argv[2]);
		cuff.roc_quant(argv[3], atof(argv[4]), atof(argv[5]));
	}

	if(string(argv[1]) == "acc-quant")
	{
		assert(argc == 5);
		gtfcuff cuff(argv[2]);
		cuff.acc_quant(argv[3], atof(argv[4]));
	}

	if(string(argv[1]) == "split")
	{
		assert(argc == 6);
		gtfcuff cuff(argv[2]);
		cuff.assign_pred(argv[3]);
		cuff.split(argv[4], argv[5]);
	}

	if(string(argv[1]) == "classify")
	{
		assert(argc == 4);
		gtfcuff cuff(argv[2]);
		cuff.assign_pred(argv[3]);
		cuff.classify();
	}

	if(string(argv[1]) == "quant")
	{
		assert(argc == 5);
		gtfcuff cuff(argv[2]);
		cuff.assign_pred(argv[3]);
		cuff.assign_ref(argv[4]);
		cuff.quant();
	}

    return 0;
}
