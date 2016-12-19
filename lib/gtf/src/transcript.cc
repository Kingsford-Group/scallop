#include <cstdio>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <map>

#include "transcript.h"

transcript::transcript()
{
}

transcript::transcript(const item &e)
{
	assign(e);
	exons.clear();
}

transcript::~transcript()
{
}

int transcript::assign(const item &e)
{
	//assert(e.feature == "transcript");
	seqname = e.seqname;
	source = e.source;
	feature = e.feature;
	gene_id = e.gene_id;
	transcript_id = e.transcript_id;
	start = e.start;
	end = e.end;
	strand = e.strand;
	frame = e.frame;
	coverage = e.coverage;
	RPKM = e.RPKM;
	FPKM = e.FPKM;
	TPM = e.TPM;
	return 0;
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

int transcript::clear()
{
	exons.clear();
	seqname = "";
	source = "";
	feature = "";
	gene_id = "";
	transcript_id = "";
	start = 0;
	end = 0;
	strand = '.';
	frame = -1;
	coverage = 0;
	RPKM = 0;
	TPM = 0;
	return 0;
}

int transcript::add_exon(int s, int t)
{
	exons.push_back(PI32(s, t));
	return 0;
}

int transcript::add_exon(const item &e)
{
	assert(e.transcript_id == transcript_id);
	add_exon(e.start, e.end);
	return 0;
}

int transcript::sort()
{
	std::sort(exons.begin(), exons.end());
	return 0;
}

int transcript::shrink()
{
	if(exons.size() == 0) return 0;
	vector<PI32> v;
	PI32 p = exons[0];
	for(int i = 1; i < exons.size(); i++)
	{
		PI32 q = exons[i];
		if(p.second == q.first)
		{
			p.second = q.second;
		}
		else
		{
			//assert(p.second < q.first);
			v.push_back(p);
			p = q;
		}
	}
	v.push_back(p);
	exons = v;
	return 0;
}

int transcript::assign_RPKM(double factor)
{
	RPKM = coverage * factor;
	return 0;
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

int transcript::write(ostream &fout) const
{
	fout.precision(4);
	fout<<fixed;

	if(exons.size() == 0) return 0;
	
	PI32 p = get_bounds();

	fout<<seqname.c_str()<<"\t";				// chromosome name
	fout<<source.c_str()<<"\t";					// source
	fout<<"transcript\t";						// feature
	fout<<p.first + 1<<"\t";					// left position
	fout<<p.second<<"\t";						// right position
	fout<<1000<<"\t";							// score, now as expression
	fout<<strand<<"\t";							// strand
	fout<<".\t";								// frame
	fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
	fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
	fout<<"RPKM \""<<RPKM<<"\"; ";
	fout<<"FPKM \""<<FPKM<<"\"; ";
	fout<<"TPM \""<<TPM<<"\"; ";
	fout<<"cov \""<<coverage<<"\";"<<endl;

	for(int k = 0; k < exons.size(); k++)
	{
		fout<<seqname.c_str()<<"\t";		// chromosome name
		fout<<source.c_str()<<"\t";			// source
		fout<<"exon\t";						// feature
		fout<<exons[k].first + 1<<"\t";		// left position
		fout<<exons[k].second<<"\t";		// right position
		fout<<1000<<"\t";					// score, now as expression
		fout<<strand<<"\t";					// strand
		fout<<".\t";						// frame
		fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
		fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
		fout<<"exon \""<<k + 1<<"\"; ";
		fout<<"RPKM \""<<RPKM<<"\"; ";
		fout<<"FPKM \""<<FPKM<<"\"; ";
		fout<<"TPM \""<<TPM<<"\"; ";
		fout<<"cov \""<<coverage<<"\";"<<endl;
	}
	return 0;
}

