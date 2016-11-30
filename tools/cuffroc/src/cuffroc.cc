#include "cuffroc.h"
#include "genome.h"
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

cuffroc::cuffroc(const string &cufffile, const string &gtffile, int r, int m)
{
	read_cuff(cufffile);
	read_gtf(gtffile);
	refsize = r;
	mexons = m;
}

int cuffroc::read_cuff(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail()) return -1;

	char line[10240];

	while(fin.getline(line, 10240, '\n'))
	{
		cuffitem x(line);
		if(x.code == '@') continue;
		items.push_back(x);
	}
	return 0;
}

int cuffroc::read_gtf(const string &file)
{
	t2e.clear();
	genome gn(file);
	for(int i = 0; i < gn.genes.size(); i++)
	{
		vector<transcript> &v = gn.genes[i].transcripts;
		for(int k = 0; k < v.size(); k++)
		{
			string s = v[k].transcript_id;
			int e = v[k].exons.size();
			if(t2e.find(s) != t2e.end()) continue;
			t2e.insert(pair<string, int>(s, e));
		}
	}
	return 0;
}

int cuffroc::filter_items()
{
	vector<cuffitem> v;
	for(int i = 0; i < items.size(); i++)
	{
		string s = items[i].transcript_id;
		if(t2e.find(s) == t2e.end()) continue;
		if(t2e[s] < mexons) continue;
		v.push_back(items[i]);
	}
	items = v;
	return 0;
}

int cuffroc::solve()
{
	filter_items();

	if(items.size() == 0) return 0;

	sort(items.begin(), items.end());

	int correct = 0;
	for(int i = 0; i < items.size(); i++) if(items[i].code == '=') correct++;

	double sen0 = correct * 100.0 / refsize;
	for(int i = 0; i < items.size(); i++)
	{
		double sen = correct * 100.0 / refsize;
		double pre = correct * 100.0 / (items.size() - i);

		if(sen * 2.0 < sen0) break;

		if(i % 100 == 0)
		{
			printf("ROC: reference = %d prediction = %lu correct = %d sensitivity = %.2lf precision = %.2lf | coverage = %.3lf\n",
				refsize, items.size() - i, correct, sen, pre, items[i].coverage);
		}

		if(items[i].code == '=') correct--;
	}

	return 0;
}

int cuffroc::print()
{
	for(int i = 0; i < items.size(); i++)
	{
		int n = 0;
		string s = items[i].transcript_id;
		if(t2e.find(s) != t2e.end()) n = t2e[s];
		items[i].print(n);
	}
	return 0;
}
