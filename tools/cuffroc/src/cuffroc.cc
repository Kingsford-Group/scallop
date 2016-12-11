#include "cuffroc.h"
#include "genome.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>

using namespace std;

cuffroc::cuffroc(const string &cufffile, const string &gtffile, int r, int m, int f, double p)
	: gm(gtffile)
{
	read_cuff(cufffile);
	build_indices();
	refsize = r;
	mexons = m;
	pratio = p;
	ftype = f;
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

	fin.close();
	return 0;
}

int cuffroc::build_indices()
{
	t2i.clear();
	for(int i = 0; i < items.size(); i++)
	{
		string s = items[i].transcript_id;
		t2i.insert(pair<string, int>(s, i));
	}

	t2e.clear();
	t2s.clear();
	for(int i = 0; i < gm.genes.size(); i++)
	{
		vector<transcript> &v = gm.genes[i].transcripts;
		for(int k = 0; k < v.size(); k++)
		{
			string s = v[k].transcript_id;
			int e = v[k].exons.size();
			char c = v[k].strand;
			if(t2e.find(s) != t2e.end()) continue;
			assert(t2s.find(s) == t2s.end());
			t2e.insert(pair<string, int>(s, e));
			t2s.insert(pair<string, char>(s, c));
		}
	}

	return 0;
}

int cuffroc::classify()
{
	ofstream f1("true.gtf");
	ofstream f2("false.gtf");
	for(int i = 0; i < gm.genes.size(); i++)
	{
		vector<transcript> &v = gm.genes[i].transcripts;
		for(int k = 0; k < v.size(); k++)
		{
			string s = v[k].transcript_id;
			bool b = false;
			if(t2i.find(s) != t2i.end())
			{
				if(items[t2i[s]].code == '=') b = true;
				else b = false;
			}
			if(b == true) v[k].write(f1);
			else v[k].write(f2);
		}
	}
	f1.close();
	f2.close();
	return 0;
}

int cuffroc::filter_items()
{
	vector<cuffitem> v;
	for(int i = 0; i < items.size(); i++)
	{
		string s = items[i].transcript_id;
		if(t2e.find(s) == t2e.end()) continue;
		if(mexons >= 0 && t2e[s] != mexons) continue;

		/* TODO
		int min_length = t2e[s] * 50 + 200;
		if(items[i].length < min_length) continue;
		if(t2e[s] == 1 && items[i].coverage < 20) continue;
		*/

		v.push_back(items[i]);
	}
	items = v;
	return 0;
}

int cuffroc::solve()
{
	filter_items();

	if(items.size() == 0) return 0;

	if(ftype == 1) sort(items.begin(), items.end(), cuffitem_cmp_coverage);
	else if(ftype == 2) sort(items.begin(), items.end(), cuffitem_cmp_length);
	else return 0;

	int correct = 0;
	for(int i = 0; i < items.size(); i++) if(items[i].code == '=') correct++;

	double max_com = 999;
	double max_sen = 0;
	double max_pre = 0;
	double max_cov = 0;
	int max_len = 0;
	int max_correct = 0;
	int max_size = 0;
	double sen0 = correct * 100.0 / refsize;
	for(int i = 0; i < items.size(); i++)
	{
		double sen = correct * 100.0 / refsize;
		double pre = correct * 100.0 / (items.size() - i);

		double com = fabs(pre - pratio);

		if(com < max_com)
		{
			max_com = com;
			max_sen = sen;
			max_pre = pre;
			max_cov = items[i].coverage;
			max_len = items[i].length;
			max_correct = correct;
			max_size = items.size() - i;
		}

		if(sen * 2.0 < sen0) break;

		if(i % 100 == 0)
		{
			printf("ROC: reference = %d prediction = %lu correct = %d sensitivity = %.2lf precision = %.2lf | coverage = %.3lf, length = %d\n",
				refsize, items.size() - i, correct, sen, pre, items[i].coverage, items[i].length);
		}

		if(items[i].code == '=') correct--;
	}

	//printf("BESTCOM: reference = %d prediction = %d correct = %d sensitivity = %.2lf precision = %.2lf | coverage = %.3lf, length = %d\n",
	//			refsize, max_size, max_correct, max_sen, max_pre, max_cov, max_len);

	return 0;
}

int cuffroc::print()
{
	for(int i = 0; i < items.size(); i++)
	{
		int n = 0;
		char c = '@';
		string s = items[i].transcript_id;
		if(t2e.find(s) != t2e.end()) n = t2e[s];
		if(t2s.find(s) != t2s.end()) c = t2s[s];
		items[i].print(n, c);
	}
	return 0;
}
