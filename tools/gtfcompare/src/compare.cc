#include "compare.h"
#include <cassert>

bool compare_transcript_structure(const transcript &x, const transcript &y)
{
	if(x.exons.size() != y.exons.size()) return false;
	for(int i = 0; i < x.exons.size(); i++)
	{
		if(x.exons[i].first != y.exons[i].first) return false;
		if(x.exons[i].second != y.exons[i].second) return false;
	}
	return true;
}

bool compare_transcript_expression(const transcript &x, const transcript &y)
{
	if(x.expression == y.expression) return true;
	else return false;
}

bool compare_transcript(const transcript &x, const transcript &y)
{
	if(compare_transcript_structure(x, y) == false) return false;
	if(compare_transcript_expression(x, y) == false) return false;
	return true;
}

int compare_gene(const gene &x, const gene &y)
{
	int cnt = 0;
	vector<bool> v;
	v.assign(y.transcripts.size(), false);
	for(int i = 0; i < x.transcripts.size(); i++)
	{
		for(int j = 0; j < y.transcripts.size(); j++)
		{
			if(v[j] == true) continue;
			if(compare_transcript(x.transcripts[i], y.transcripts[j]) == false) continue;
			v[j] = true;
			cnt++;
			break;
		}
	}
	return cnt;
}

int compare_genome(const genome &x, const genome &y)
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
		int t0 = compare_gene(*gx, *gy);
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
