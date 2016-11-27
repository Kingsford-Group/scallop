#include "compare.h"
#include "config.h"
#include <cassert>
#include <cmath>
#include <algorithm>

bool compare_structure(const transcript &x, const transcript &y)
{
	if(x.exons.size() != y.exons.size()) return false;
	for(int i = 0; i < x.exons.size(); i++)
	{
		if(x.exons[i].first != y.exons[i].first) return false;
		if(x.exons[i].second != y.exons[i].second) return false;
	}
	return true;
}

bool compare_intron_chain(const transcript &x, const transcript &y)
{
	if(x.exons.size() != y.exons.size()) return false;
	
	if(x.exons.size() == 0) return false;

	if(x.exons.size() == 1)
	{
		if(fabs(x.exons[0].first - y.exons[0].first) > 100) return false;
		if(fabs(x.exons[0].second - y.exons[0].second) > 100) return false;
		return true;
	}

	for(int i = 0; i < x.exons.size(); i++)
	{
		double diff1 = 0.5;
		double diff2 = 0.5;
		if(i == 0) diff1 = 9999999999;
		if(i == x.exons.size() - 1) diff2 = 9999999999;

		if(fabs(x.exons[i].first - y.exons[i].first) > diff1) return false;
		if(fabs(x.exons[i].second - y.exons[i].second) > diff2) return false;
	}
	return true;
}

bool compare_expression(const transcript &x, const transcript &y)
{
	if(fabs(x.coverage - y.coverage) < 0.5) return true;
	else return false;
}

int compare_transcripts(const vector<transcript> &y, const vector<transcript> &x)
{
	vector<bool> v;
	return compare_transcripts(y, x, v);
}

int compare_transcripts(const vector<transcript> &y, const vector<transcript> &x, vector<bool> &vv)
{
	int cnt = 0;
	vector<bool> v;
	v.assign(y.size(), false);
	vv.clear();
	for(int i = 0; i < x.size(); i++)
	{
		const transcript &t1 = x[i];
		bool flag = false;
		for(int j = 0; j < y.size(); j++)
		{
			if(v[j] == true) continue;
			const transcript &t2 = y[j];
			bool b = false;
			if(algo == 1)
			{
				if(compare_structure(t1, t2) == false) b = false;
				if(compare_expression(t1, t2) == false) b = false;
			}
			if(algo == 2)
			{
				b = compare_intron_chain(t1, t2);
				if(t1.strand != t2.strand) b = false;
			}

			if(b == false) continue;

			printf("TRUE %s %s %s %c %lu %d %.2lf %.3lf %.2lf %s %s %s %c %lu %d %.2lf %.3lf %.2lf\n", 
					t1.gene_id.c_str(), t1.transcript_id.c_str(), t1.label().c_str(), t1.strand, 
					t1.exons.size(), t1.length(), t1.coverage, t1.covratio, t1.RPKM,
					t2.gene_id.c_str(), t2.transcript_id.c_str(), t2.label().c_str(), t2.strand, 
					t2.exons.size(), t2.length(), t2.coverage, t2.covratio, t2.RPKM);
			
			flag = true;
			v[j] = true;
			cnt++;
			break;
		}

		vv.push_back(flag);

		if(flag == false)
		{
			printf("FALSE %s %s %s %c %lu %d %.2lf %.3lf %.2lf\n", 
				t1.gene_id.c_str(), t1.transcript_id.c_str(), t1.label().c_str(), t1.strand, 
				t1.exons.size(), t1.length(), t1.coverage, t1.covratio, t1.RPKM);

		}
	}
	return cnt;
}

int compare_genome1(const genome &x, const genome &y)
{
	int gequal = 0;
	int tequal = 0;
	int gtotal = x.genes.size();
	int ttotal = 0;
	for(int i = 0; i < x.genes.size(); i++)
	{
		const gene* gx = &(x.genes[i]);
		const gene* gy = y.get_gene(x.genes[i].get_gene_id());
		if(gx == NULL || gy == NULL) continue;
		int tx = gx->transcripts.size();
		int ty = gy->transcripts.size();
		int t0 = compare_transcripts(gx->transcripts, gy->transcripts);
		assert(t0 <= tx);
		assert(t0 <= ty);
		ttotal += tx;
		tequal += t0;
		if(t0 == tx) gequal++;
		string s;
		if(tx == ty) s = "EQUAL";
		else if(tx > ty) s = "GREATER";
		else s = "LESS";
		printf("%s %d and %d transcripts, %d are equal, %s, %s\n", gx->get_gene_id().c_str(), tx, ty, t0, (t0 == tx) ? "TRUE" : "FALSE", s.c_str());
	}
	printf("summary: %d out of %d genes are equal, %d out of %d transcripts are equal\n",
			gequal, gtotal, tequal, ttotal);
	return 0;
}

int compare_genome2(const genome &x, const genome &y)
{
	typedef pair< string, vector<transcript> > PSVT;
	typedef map< string, vector<transcript> > MSVT;
	MSVT m1;
	MSVT m2;

	int xtotal = 0;
	int ytotal = 0;
	for(int i = 0; i < x.genes.size(); i++)
	{
		string chrm = x.genes[i].get_seqname();
		const vector<transcript> &v0 = x.genes[i].transcripts;
		vector<transcript> v;
		for(int k = 0; k < v0.size(); k++)
		{
			if(v0[k].exons.size() <= 1 && multiple_exon == true) continue;
			v.push_back(v0[k]);
		}
		xtotal += v.size();
		if(m1.find(chrm) == m1.end())
		{
			m1.insert(PSVT(chrm, v));
		}
		else
		{
			m1[chrm].insert(m1[chrm].end(), v.begin(), v.end());
		}
	}

	for(int i = 0; i < y.genes.size(); i++)
	{
		string chrm = y.genes[i].get_seqname();
		const vector<transcript> &v0 = y.genes[i].transcripts;
		vector<transcript> v;
		for(int k = 0; k < v0.size(); k++)
		{
			if(v0[k].exons.size() <= 1 && multiple_exon == true) continue;
			if(v0[k].length() < min_transcript_length) continue;
			v.push_back(v0[k]);
		}
		ytotal += v.size();
		if(m2.find(chrm) == m2.end())
		{
			m2.insert(PSVT(chrm, v));
		}
		else
		{
			m2[chrm].insert(m2[chrm].end(), v.begin(), v.end());
		}
	}

	int correct = 0;
	vector<PTB> vv;
	for(MSVT::iterator it = m1.begin(); it != m1.end(); it++)
	{
		const vector<transcript> &v1 = it->second;
		if(m2.find(it->first) == m2.end()) continue;
		const vector<transcript> &v2 = m2[it->first];
		vector<bool> vb;
		correct += compare_transcripts(v1, v2, vb);
		assert(v2.size() == vb.size());
		for(int k = 0; k < v2.size(); k++) vv.push_back(PTB(v2[k], vb[k]));
	}

	double s = (xtotal == 0) ? 0 : correct * 100.0 / xtotal;
	double p = (ytotal == 0) ? 0 : correct * 100.0 / ytotal;
	printf("SUMMARY: reference = %d prediction = %d correct = %d sensitivity = %.2lf precision = %.2lf\n", xtotal, ytotal, correct, s, p);

	ROC(vv, xtotal);

	return 0;
}

int ROC(vector<PTB> &vv, int xtotal)
{
	if(vv.size() == 0) return 0;

	sort(vv.begin(), vv.end(), transcript_cmp);
	int correct = 0;
	for(int i = 0; i < vv.size(); i++) if(vv[i].second == true) correct++;

	double sen0 = correct * 100.0 / xtotal;
	for(int i = 0; i < vv.size(); i++)
	{
		double sen = correct * 100.0 / xtotal;
		double pre = correct * 100.0 / (vv.size() - i);

		if(sen * 2.0 < sen0) break;

		if(i % 100 == 0)
		{
			printf("ROC: reference = %d prediction = %lu correct = %d sensitivity = %.2lf precision = %.2lf | coverage = %.3lf\n",
				xtotal, vv.size() - i, correct, sen, pre, vv[i].first.coverage);
		}

		if(vv[i].second == true) correct--;
	}
	return 0;
}

bool transcript_cmp(const PTB &x, const PTB &y)
{
	if(x.first.coverage < y.first.coverage) return true;
	else return false;
}
