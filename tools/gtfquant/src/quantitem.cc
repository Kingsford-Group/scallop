#include "quantitem.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <cassert>
#include <cstdio>
#include <cmath>

quantitem::quantitem(const string &s)
{
	parse(s);
}

int quantitem::parse(const string &s)
{
	char buf[10240];
	stringstream sstr(s);
	sstr>>buf;
	transcript_id.assign(buf);
	sstr>>length;
	sstr>>elength;
	sstr>>tpm;
	sstr>>numreads;
	return 0;
}
