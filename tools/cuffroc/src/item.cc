#include "item.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <cassert>
#include <cstdio>
#include <cmath>

item::item(const string &s)
{
	code = '@';
	parse(s);
}

int item::parse(const string &s)
{
	char buf[10240];
	stringstream sstr(s);
	sstr>>buf;
	if(string(buf) == "ref_gene_id") return 0;
	sstr>>buf>>buf;
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

int item::print() const
{
	printf("%s %s %c %.3lf %d\n", transcript_id.c_str(), gene_id.c_str(), code, coverage, length);
	return 0;
}

bool item::operator<(const item &ge) const
{
	if(coverage < ge.coverage) return true;
	else return false;
}

