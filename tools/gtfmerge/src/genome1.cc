#include "genome1.h"
#include "config.h"
#include <cassert>
#include <fstream>

genome1::genome1()
{}

genome1::genome1(const string &file)
{
	build(file);
}

int genome1::build(const string &file)
{
	build_multiexon_transcripts(file);
	build_intron_index();
	return 0;
}

int genome1::build(const vector<transcript> &v)
{
	transcripts = v;
	build_intron_index();
	return 0;
}

int genome1::clear()
{
	transcripts.clear();
	intron_index.clear();
	return 0;
}

int genome1::build_multiexon_transcripts(const string &file)
{
	genome gm(file);
	for(int i = 0; i < gm.genes.size(); i++)
	{
		const gene &g = gm.genes[i];
		for(int k = 0; k < g.transcripts.size(); k++)
		{
			const transcript &t = g.transcripts[k];
			if(t.exons.size() <= 1) continue;
			transcripts.push_back(t);
		}
	}
	return 0;
}

int genome1::build_intron_index()
{
	intron_index.clear();
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcript &t = transcripts[i];
		assert(t.exons.size() >= 2);
		PI32 p = t.get_first_intron();
		if(intron_index.find(p) == intron_index.end()) 
		{
			set<int> s;
			s.insert(i);
			intron_index.insert(PPIS(p, s));
		}
		else
		{
			intron_index[p].insert(i);
		}
	}
	return 0;
}

int genome1::query(const transcript &t, const set<int> &fb)
{
	if(t.exons.size() <= 1) return -1;
	PI32 p = t.get_first_intron();
	if(intron_index.find(p) == intron_index.end()) return -1;
	set<int> s = intron_index[p];
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		transcript &x = transcripts[k];
		if(x.strand != t.strand) continue;
		if(x.seqname != t.seqname) continue;
		if(x.exons.size() != t.exons.size()) continue;
		if(x.intron_chain_match(t) == false) continue;
		if(fb.find(k) != fb.end()) continue;
		return k;
	}
	return -1;
}

int genome1::compare(const genome1 &gy, MII &x2y, MII &y2x)
{
	x2y.clear();
	y2x.clear();
	set<int> fb;
	for(int i = 0; i < gy.transcripts.size(); i++)
	{
		const transcript &t = gy.transcripts[i];
		int k = query(t, fb);
		if(k == -1) continue;
		x2y.insert(PII(k, i));
		y2x.insert(PII(i, k));
		fb.insert(k);
	}
	return 0;
}

int genome1::build_intersection(const genome1 &gm, genome1 &out)
{
	vector<transcript> vv;
	MII x2y;
	MII y2x;
	compare(gm, x2y, y2x);
	for(MII::iterator it = x2y.begin(); it != x2y.end(); it++)
	{
		int i = it->first;
		int j = it->second;
		transcript x = transcripts[i];
		const transcript &y = gm.transcripts[j];
		x.coverage += y.coverage;
		x.RPKM += y.RPKM;
		x.FPKM += y.FPKM;
		x.TPM += y.TPM;
		vv.push_back(x);
	}

	out.clear();
	out.build(vv);
	return 0;
}

int genome1::build_union(const genome1 &gm)
{
	MII x2y;
	MII y2x;
	compare(gm, x2y, y2x);
	for(MII::iterator it = x2y.begin(); it != x2y.end(); it++)
	{
		int i = it->first;
		int j = it->second;
		transcript &x = transcripts[i];
		const transcript &y = gm.transcripts[j];
		x.coverage += y.coverage;
		x.RPKM += y.RPKM;
		x.FPKM += y.FPKM;
		x.TPM += y.TPM;
	}

	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		if(y2x.find(i) != y2x.end()) continue;
		transcripts.push_back(gm.transcripts[i]);
	}

	build_intron_index();
	return 0;
}

int genome1::add_suffix(const string &p)
{
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcript &t = transcripts[i];
		t.transcript_id.append("-").append(p);
		t.gene_id.append("-").append(p);
	}
	return 0;
}

int genome1::print(int index)
{
	printf("genome %d: %lu transcripts, %lu distinct first intron\n", index, transcripts.size(), intron_index.size());
	return 0;
}

int genome1::write(const string &file)
{
	ofstream fout(file.c_str());
	if(fout.fail()) return 0;
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcript &t = transcripts[i];
		if(t.coverage < min_transcript_coverage) continue;
		t.write(fout);
	}
	fout.close();
}
