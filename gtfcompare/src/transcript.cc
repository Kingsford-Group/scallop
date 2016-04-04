#include <cstdio>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <map>

#include "transcript.h"

transcript::transcript()
{
}

transcript::~transcript()
{
}

int transcript::add_exon(const exon &e)
{
	seqname = e.seqname;
	source = e.source;
	feature = e.feature;
	transcript_id = e.transcript_id;
	gene_id = e.gene_id;
	expression = e.expression;
	exons.push_back(PI32(e.start, e.end));
	return 0;
}

int transcript::sort()
{
	std::sort(exons.begin(), exons.end());
	return 0;
}

bool transcript::equal_structure(const transcript &t) const
{
	if(t.exons.size() != exons.size()) return false;
	for(int i = 0; i < exons.size(); i++)
	{
		if(exons[i].first != t.exons[i].first) return false;
		if(exons[i].second != t.exons[i].second) return false;
	}
	return true;
}

bool transcript::equal_expression(const transcript &t) const
{
	if(t.expression == expression) return true;
	else return false;
}

bool transcript::equal(const transcript &t) const
{
	if(equal_structure(t) == false) return false;
	if(equal_expression(t) == false) return false;
	return true;
}

int transcript::write(ofstream &fout) const
{
	fout.precision(2);
	fout<<fixed;

	if(exons.size() == 0) return 0;

	int lpos = exons[0].first;
	int rpos = exons[exons.size() - 1].second;

	fout<<seqname.c_str()<<"\t";				// chromosome name
	fout<<source.c_str()<<"\t";					// source
	fout<<"transcript\t";						// feature
	fout<<lpos<<"\t";							// left position
	fout<<rpos<<"\t";							// right position
	fout<<1000<<"\t";							// score, now as expression
	fout<<"+\t";								// strand
	fout<<".\t";								// frame
	fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
	fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
	fout<<"expression \""<<expression<<"\";"<<endl;

	for(int k = 0; k < exons.size(); k++)
	{
		fout<<seqname.c_str()<<"\t";		// chromosome name
		fout<<source.c_str()<<"\t";			// source
		fout<<"exon\t";						// feature
		fout<<exons[k].first<<"\t";			// left position
		fout<<exons[k].second - 1<<"\t";	// right position
		fout<<1000<<"\t";					// score, now as expression
		fout<<"+\t";						// strand
		fout<<".\t";						// frame
		fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
		fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
		fout<<"exon \""<<k + 1<<"\"; ";
		fout<<"expression \""<<expression<<"\";"<<endl;
	}
	return 0;
}

