#include "exon.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <cassert>
#include <cstdio>

exon::exon(const string &s)
{
	parse(s);
}

exon::exon(const string &_transcript_id, const string &_gene_id, int32_t _start, int32_t _end, int32_t _expression)
	: seqname("chr1"), source("simulation"), start(_start), end(_end), transcript_id(_transcript_id),
	gene_id(_gene_id), expression(_expression), feature("exon"), score(1000), strand('+'), frame(0)
{
}


int exon::parse(const string &s)
{
	char buf[10240];
	stringstream sstr(s);
	sstr>>buf;
	seqname.assign(buf);
	sstr>>buf;
	source.assign(buf);
	sstr>>buf;
	feature.assign(buf);
	sstr>>start>>end;
	end++;			// TODO check gtf end is inclusive
	sstr>>buf;
	if(buf[0] == '.') score = -1;
	else score = atof(buf);
	sstr>>buf;
	strand = buf[0];
	sstr>>buf;
	frame = buf[0];

	char buf2[10240];
	while(sstr.eof() == false)
	{
		sstr>>buf;
		sstr.getline(buf2, 10240, ';');
		string k(buf2);
		if(string(buf) == "" || k == "") break;

		size_t p1 = k.find_first_of('"');
		size_t p2 = k.find_last_of('"');
		assert(p1 != string::npos);
		assert(p2 != string::npos);
		assert(p1 != p2);
		string v = k.substr(p1 + 1, p2 - p1 - 1);

		//printf(" |%s|%s|\n", buf, v.c_str());

		if(string(buf) == "transcript_id") transcript_id = v;
		else if(string(buf) == "gene_id") gene_id = v;
		else if(string(buf) == "expression") expression = atoi(v.c_str());
	}
	
	return 0;
}

int exon::print() const
{
	printf("%s\t%s\t%s\t%d\t%d\t%.1lf\t%c\t%c\ttranscript_id \"%s\"; gene_id \"%s\"; expression \"%d\"\n",
			seqname.c_str(), source.c_str(), feature.c_str(), start, end, score, strand, frame,
			transcript_id.c_str(), gene_id.c_str(), expression);
	return 0;
}

bool exon::operator<(const exon &ge) const
{
	if(start < ge.start) return true;
	else return false;
}

int exon::length() const
{
	assert(end > start);
	return end - start;
}
