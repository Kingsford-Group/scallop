#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include "genome.h"

using namespace std;

int main(int argc, const char **argv)
{
 	if(argc == 1)
	{
		cout<<"usage: " << endl;
		cout<<"       " << argv[0] << " RPKM2TPM <in-gtf-file> <out-gtf-file>"<<endl;
		cout<<"       " << argv[0] << " FPKM2TPM <in-gtf-file> <out-gtf-file>"<<endl;
		cout<<"       " << argv[0] << " format <in-gtf-file> <out-gtf-file>"<<endl;
		cout<<"       " << argv[0] << " filter <min-transcript-coverage> <in-gtf-file> <out-gtf-file>"<<endl;
		return 0;
	}

	if(string(argv[1]) == "FPKM2TPM")
	{
		genome gm(argv[2]);
		gm.assign_TPM_by_FPKM();
		gm.write(argv[3]);
	}

	if(string(argv[1]) == "RPKM2TPM")
	{
		genome gm(argv[2]);
		gm.assign_TPM_by_RPKM();
		gm.write(argv[3]);
	}

	if(string(argv[1]) == "format")
	{
		genome gm(argv[2]);
		gm.write(argv[3]);
	}

	if(string(argv[1]) == "filter")
	{
		double c = atof(argv[2]);
		genome gm(argv[3]);
		gm.filter_low_coverage_transcripts(c);
		gm.write(argv[4]);
	}

    return 0;
}
