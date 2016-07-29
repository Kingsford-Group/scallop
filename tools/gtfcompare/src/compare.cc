#include "compare.h"
#include <cassert>

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
	for(int i = 0; i < x.exons.size(); i++)
	{
		if(i >= 1 && x.exons[i].first != y.exons[i].first) return false;
		if(i < x.exons.size() - 1 && x.exons[i].second != y.exons[i].second) return false;
	}
	return true;
}

bool compare_expression(const transcript &x, const transcript &y)
{
	if(x.expression == y.expression) return true;
	else return false;
}

int compare_gene(const gene &x, const gene &y, int mode)
{
	return compare_transcripts(x.transcripts, y.transcripts, mode);
}

int compare_transcripts(const vector<transcript> &x, const vector<transcript> &y, int mode)
{
	int cnt = 0;
	vector<bool> v;
	v.assign(y.size(), false);
	for(int i = 0; i < x.size(); i++)
	{
		const transcript &t1 = x[i];
		for(int j = 0; j < y.size(); j++)
		{
			if(v[j] == true) continue;
			const transcript &t2 = y[j];
			bool b = false;
			if(mode == 1)
			{
				if(compare_structure(t1, t2) == false) b = false;
				if(compare_expression(t1, t2) == false) b = false;
			}
			if(mode == 2)
			{
				b = compare_intron_chain(t1, t2);
				if(t1.strand != t2.strand) b = false;
			}

			if(b == false) continue;

			v[j] = true;
			cnt++;
			break;
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
		int t0 = compare_gene(*gx, *gy, 1);
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
		const vector<transcript> &v = x.genes[i].transcripts;
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
		const vector<transcript> &v = y.genes[i].transcripts;
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
	for(MSVT::iterator it = m1.begin(); it != m1.end(); it++)
	{
		const vector<transcript> &v1 = it->second;
		if(m2.find(it->first) == m2.end()) continue;
		const vector<transcript> &v2 = m2[it->first];

		correct += compare_transcripts(v1, v2, 2);
	}

	double s = (xtotal == 0) ? 0 : correct * 100.0 / xtotal;
	double p = (ytotal == 0) ? 0 : correct * 100.0 / ytotal;
	printf("reference = %d prediction = %d correct = %d sensitivity = %.2lf precision = %.2lf\n", xtotal, ytotal, correct, s, p);

	return 0;
}
