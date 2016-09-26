#include <cstdio>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <map>

#include "transcript.h"

transcript::transcript()
{
	RPKM = 0;
}

transcript::~transcript()
{
}

bool transcript::operator< (const transcript &t) const
{
	int b = seqname.compare(t.seqname);
	if(b < 0) return true;
	if(b > 0) return false;
	if(exons.size() == 0) return true;
	if(t.exons.size() == 0) return false;
	if(exons[0].first < t.exons[0].first) return true;
	else return false;
}

int transcript::add_exon(int s, int t)
{
	exons.push_back(PI32(s, t));
	return 0;
}

int transcript::add_exon(const exon &e)
{
	seqname = e.seqname;
	source = e.source;
	feature = e.feature;
	transcript_id = e.transcript_id;
	gene_id = e.gene_id;
	expression = e.expression;
	coverage = e.coverage;
	strand = e.strand;
	exons.push_back(PI32(e.start, e.end));
	return 0;
}

int transcript::sort()
{
	std::sort(exons.begin(), exons.end());
	return 0;
}

int transcript::assign_RPKM(int reads)
{
	RPKM = coverage * 1e3 / length() * 1e6 / reads;
}

int transcript::length() const
{
	int s = 0;
	for(int i = 0; i < exons.size(); i++)
	{
		assert(exons[i].second > exons[i].first);
		s += exons[i].second - exons[i].first;
	}
	return s;
}

PI32 transcript::get_bounds() const
{
	if(exons.size() == 0) return PI32(-1, -1);
	int32_t p = exons[0].first;
	int32_t q = exons[exons.size() - 1].second;
	return PI32(p, q);
}

string transcript::label() const
{
	char buf[10240];
	PI32 p = get_bounds();
	sprintf(buf, "%s:%d-%d", seqname.c_str(), p.first, p.second);
	return string(buf);
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
	fout<<lpos + 1<<"\t";							// left position
	fout<<rpos<<"\t";							// right position
	fout<<1000<<"\t";							// score, now as expression
	fout<<strand<<"\t";					// strand
	fout<<".\t";								// frame
	fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
	fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
	fout<<"RPKM \""<<RPKM<<"\"; ";
	fout<<"coverage \""<<coverage<<"\"; ";
	fout<<"expression \""<<expression<<"\";"<<endl;

	for(int k = 0; k < exons.size(); k++)
	{
		fout<<seqname.c_str()<<"\t";		// chromosome name
		fout<<source.c_str()<<"\t";			// source
		fout<<"exon\t";						// feature
		fout<<exons[k].first + 1<<"\t";			// left position
		fout<<exons[k].second<<"\t";	// right position
		fout<<1000<<"\t";					// score, now as expression
		fout<<strand<<"\t";					// strand
		fout<<".\t";						// frame
		fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
		fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
		fout<<"exon \""<<k + 1<<"\"; ";
		fout<<"RPKM \""<<RPKM<<"\"; ";
		fout<<"coverage \""<<coverage<<"\"; ";
		fout<<"expression \""<<expression<<"\";"<<endl;
	}
	return 0;
}

