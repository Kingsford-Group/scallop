#include "cuffitem.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <cassert>
#include <cstdio>
#include <cmath>

cuffitem::cuffitem(const string &s)
{
	code = '@';
	parse(s);
}

int cuffitem::parse(const string &s)
{
	char buf[10240];
	stringstream sstr(s);
	sstr>>buf;
	ref_gene_id.assign(buf);
	if(ref_gene_id == "ref_gene_id") return 0;
	sstr>>buf;
	ref_transcript_id.assign(buf);
	sstr>>buf;
	code = buf[0];
	sstr>>buf;
	gene_id.assign(buf);
	sstr>>buf;
	transcript_id.assign(buf);
	sstr>>buf>>buf>>buf>>buf>>buf;
	coverage = atof(buf);
	sstr>>buf;
	length = atoi(buf);
	return 0;
}

int cuffitem::print(int n) const
{
	printf("%s %s %s %s %c %.3lf %d %d\n", ref_gene_id.c_str(), ref_transcript_id.c_str(), transcript_id.c_str(), gene_id.c_str(), code, coverage, length, n);
	return 0;
}

bool cuffitem::operator<(const cuffitem &ge) const
{
	if(coverage < ge.coverage) return true;
	//if(length < ge.length) return true;
	else return false;
}

