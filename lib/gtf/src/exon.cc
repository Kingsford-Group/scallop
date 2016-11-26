#include "exon.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <cassert>
#include <cstdio>
#include <cmath>

exon::exon(const string &s)
{
	parse(s);
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
	start--;			// gtf: (from 1, both inclusive)
	sstr>>buf;
	if(buf[0] == '.') score = -1;
	else score = atof(buf);
	sstr>>buf;
	strand = buf[0];
	sstr>>buf;
	frame = buf[0];

	char buf2[10240];
	expression = 0;
	coverage = 0;
	while(sstr.eof() == false)
	{
		sstr>>buf;
		/*
		sstr.getline(buf2, 10240, '"');
		sstr.getline(buf2, 10240, '"');
		*/
		sstr.getline(buf2, 10240, ';');
		string v(buf2);
		int k1 = v.find('"');
		int k2 = v.rfind('"');
		if(k1 >= 0 && k1 < k2 && k2 < v.size()) v = v.substr(k1 + 1, k2 - k1 - 1);

		if(string(buf) == "" || v == "") break;

		//printf("|%s|%s|\n", buf, v.c_str());

		if(string(buf) == "transcript_id") transcript_id = v;
		else if(string(buf) == "gene_id") gene_id = v;
		else if(string(buf) == "expression") expression = atoi(v.c_str());
		else if(string(buf) == "cov") coverage = atof(v.c_str());
		else if(string(buf) == "coverage") coverage = atof(v.c_str());
		else if(string(buf) == "numreads") numreads = atoi(v.c_str());
		else if(string(buf) == "RPKM") RPKM = atof(v.c_str());

		//sstr.getline(buf2, 10240, ';');
	}

	return 0;
}

int exon::print() const
{
	printf("%s\t%s\t%s\t%d\t%d\t%.1lf\t%c\t%c\ttranscript_id \"%s\"; gene_id \"%s\"; expression \"%d\"; coverage \"%.2f\"; numreads \"%d\"; RPKM \"%.2lf\"\n",
			seqname.c_str(), source.c_str(), feature.c_str(), start, end, score, strand, frame,
			transcript_id.c_str(), gene_id.c_str(), expression, coverage, numreads, RPKM);
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

bool exon::merge(const exon &e)
{
	if(end != e.start) return false;
	assert(seqname == e.seqname);
	assert(source == e.source);
	assert(feature == e.feature);
	assert(transcript_id == e.transcript_id);
	assert(gene_id == e.gene_id);
	assert(fabs(coverage - e.coverage) < 0.1);
	assert(fabs(RPKM - e.RPKM) < 0.1);
	assert(numreads == numreads);
	assert(expression == e.expression);
	assert(frame == e.frame);
	assert(strand == e.strand);
	end = e.end;
	return true;
}
