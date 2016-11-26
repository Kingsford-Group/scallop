#include <cstdlib>
#include <algorithm>
#include "simulator.h"
#include "util.h"

simulator::simulator(int _num_exons, int _num_transcripts, int _max_length, int _max_expression)
	: num_exons(_num_exons), num_transcripts(_num_transcripts), max_expression(_max_expression), max_length(_max_length)
{
}

int simulator::simulate_transcript(const string &tid, const string &gid, transcript &t)
{
	// number of exons in this transcript
	int num = rand() % max_length + 1;

	// gene num numbers from [0, ne)
	vector<int> v;
	for(int i = 0; i < num_exons; i++) v.push_back(i);

	vector<int> s;
	for(int i = 0; i < num; i++)
	{
		int r = rand() % (num_exons - i);
		s.push_back(v[r]);
		int p = num_exons - i - 1;
		int t = v[r];
		v[r] = v[p];
		v[p] = t;
	}
	sort(s.begin(), s.end());

	// generate expression level
	int x = rand() % max_expression + 1;

	t.clear();
	t.transcript_id = tid;
	t.gene_id = gid;
	t.coverage = x;

	// build transcript
	for(int i = 0; i < s.size(); i++)
	{
		int ss = s[i] * 3 + 1;
		int tt = s[i] * 3 + 3;
		t.add_exon(ss, tt);
	}
	return 0;
}

int simulator::simulate_gene(const string &gid, gene &g)
{
	for(int i = 0; i < num_transcripts; i++)
	{
		string tid = "transcript" + tostring(i + 1);
		transcript t;
		simulate_transcript(tid, gid, t);
		g.add_transcript(t);
	}
	return 0;
}
