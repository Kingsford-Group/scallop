#include "compare.h"
#include <cassert>
#include <cmath>

int remove_single_exon_transcripts(genome &gm)
{
	for(int i = 0; i < gm.genes.size(); i++)
	{
		gm.genes[i].remove_single_exon_transcripts();
	}
	return 0;
}

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
	if(x.expression == y.expression) return true;
	else return false;
}

int compare_gene(const gene &x, const gene &y, int mode)
{
	return compare_transcripts(x.transcripts, y.transcripts, mode);
}

int compare_transcripts(const vector<transcript> &y, const vector<transcript> &x, int mode)
{
	int cnt = 0;
	vector<bool> v;
	v.assign(y.size(), false);
	for(int i = 0; i < x.size(); i++)
	{
		const transcript &t1 = x[i];
		bool flag = false;
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

				if(b == true)
				{
					printf("TRUE %s %s %s %c %s %s %s %c\n", t1.gene_id.c_str(), t1.transcript_id.c_str(), t1.label().c_str(), t1.strand, 
							t2.gene_id.c_str(), t2.transcript_id.c_str(), t2.label().c_str(), t2.strand);
				}

				if(t1.strand != t2.strand) b = false;
			}

			if(b == false) continue;

			flag = true;
			v[j] = true;
			cnt++;
			break;
		}

		if(flag == false) printf("FALSE %s %s %s %c\n", t1.gene_id.c_str(), t1.transcript_id.c_str(), t1.label().c_str(), t1.strand);
	}
	return cnt;
}

int compare_gene_bounds(const gene &x, const gene &y)
{
	if(x.get_seqname() != y.get_seqname()) return 0;
	PI32 px = x.get_bounds();
	PI32 py = y.get_bounds();
	assert(px.first < px.second);
	assert(py.first < py.second);

	if(py.first < px.first && px.first < py.second && py.second < px.second) printf("%s %c %s:%d-%d %s %c %s:%d-%d overlap1\n", x.get_gene_id().c_str(), x.get_strand(), x.get_seqname().c_str(), px.first, px.second, y.get_gene_id().c_str(), y.get_strand(), y.get_seqname().c_str(), py.first, py.second);
	if(px.first < py.first && py.first < px.second && px.second < py.second) printf("%s %c %s:%d-%d %s %c %s:%d-%d overlap2\n", x.get_gene_id().c_str(), x.get_strand(), x.get_seqname().c_str(), px.first, px.second, y.get_gene_id().c_str(), y.get_strand(), y.get_seqname().c_str(), py.first, py.second);
	if(px.first <= py.first && px.second >= py.second) printf("%s %c %s:%d-%d %s %c %s:%d-%d inclusive1\n", x.get_gene_id().c_str(), x.get_strand(), x.get_seqname().c_str(), px.first, px.second, y.get_gene_id().c_str(), y.get_strand(), y.get_seqname().c_str(), py.first, py.second);
	if(py.first <= px.first && py.second >= px.second) printf("%s %c %s:%d-%d %s %c %s:%d-%d inclusive2\n", x.get_gene_id().c_str(), x.get_strand(), x.get_seqname().c_str(), px.first, px.second, y.get_gene_id().c_str(), y.get_strand(), y.get_seqname().c_str(), py.first, py.second);

	return 0;
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

int compare_genome3(const genome &x, const genome &y)
{
	typedef pair< string, vector<int> > PSVI;
	typedef map< string, vector<int> > MSVI;
	MSVI m1;
	MSVI m2;

	for(int i = 0; i < x.genes.size(); i++)
	{
		string chrm = x.genes[i].get_seqname();
		if(m1.find(chrm) == m1.end())
		{
			vector<int> v;
			v.push_back(i);
			m1.insert(PSVI(chrm, v));
		}
		else
		{
			m1[chrm].push_back(i);
		}
	}

	for(int i = 0; i < y.genes.size(); i++)
	{
		string chrm = y.genes[i].get_seqname();
		if(m2.find(chrm) == m2.end())
		{
			vector<int> v;
			v.push_back(i);
			m2.insert(PSVI(chrm, v));
		}
		else
		{
			m2[chrm].push_back(i);
		}
	}

	for(MSVI::iterator it = m1.begin(); it != m1.end(); it++)
	{
		vector<int> v1 = it->second;
		if(m2.find(it->first) == m2.end()) continue;
		vector<int> v2 = m2[it->first];

		for(int i = 0; i < v1.size(); i++)
		{
			for(int j = 0; j < v2.size(); j++)
			{
				compare_gene_bounds(x.genes[v1[i]], y.genes[v2[j]]);
			}
		}
	}

	return 0;
}

