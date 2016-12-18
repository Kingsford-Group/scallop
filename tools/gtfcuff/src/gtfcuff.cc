#include "gtfcuff.h"
#include "genome.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>

using namespace std;

gtfcuff::gtfcuff(const string &cufffile, const string &gtffile)
	: gm(gtffile)
{
	read_cuff(cufffile);
	build_indices();
	refsize = 0;
}

int gtfcuff::read_cuff(const string &file)
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

int gtfcuff::build_indices()
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

int gtfcuff::classify(const string &fn1, const string &fn2)
{
	ofstream f1(fn1.c_str());
	ofstream f2(fn2.c_str());
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

int gtfcuff::filter_items()
{
	vector<cuffitem> v;
	for(int i = 0; i < items.size(); i++)
	{
		string s = items[i].transcript_id;
		if(t2e.find(s) == t2e.end()) continue;

		/* TODO
		if(mexons >= 0 && t2e[s] != mexons) continue;
		int min_length = t2e[s] * 50 + 200;
		if(items[i].length < min_length) continue;
		if(t2e[s] == 1 && items[i].coverage < 20) continue;
		*/

		v.push_back(items[i]);
	}
	items = v;
	return 0;
}

int gtfcuff::roc()
{
	if(items.size() == 0) return 0;

	/*
	if(ftype == 1) sort(items.begin(), items.end(), cuffitem_cmp_coverage);
	else if(ftype == 2) sort(items.begin(), items.end(), cuffitem_cmp_length);
	else return 0;
	*/

	int correct = 0;
	for(int i = 0; i < items.size(); i++) if(items[i].code == '=') correct++;

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

		if(sen * 2.0 < sen0) break;

		if(i % 100 == 0)
		{
			printf("ROC: reference = %d prediction = %lu correct = %d sensitivity = %.2lf precision = %.2lf | coverage = %.3lf, length = %d\n",
				refsize, items.size() - i, correct, sen, pre, items[i].coverage, items[i].length);
		}

		if(items[i].code == '=') correct--;
	}

	return 0;
}

int gtfcuff::print()
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
