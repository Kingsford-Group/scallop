#include "gtfquant.h"
#include "genome.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>

using namespace std;

gtfquant::gtfquant(const string &quantfile, const string &gtffile, double m1, double m2)
	: gm(gtffile)
{
	read_quant(quantfile);
	build_indices();
	min_tpm = m1;
	min_numreads = m2;
}

int gtfquant::read_quant(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail()) return -1;

	char line[10240];
	while(fin.getline(line, 10240, '\n'))
	{
		quantitem x(line);
		items.push_back(x);
	}

	fin.close();
	return 0;
}

int gtfquant::build_indices()
{
	t2i.clear();
	for(int i = 0; i < items.size(); i++)
	{
		string s = items[i].transcript_id;
		t2i.insert(pair<string, int>(s, i));
	}
	return 0;
}

int gtfquant::filter()
{
	for(int i = 0; i < gm.genes.size(); i++)
	{
		vector<transcript> &v = gm.genes[i].transcripts;
		for(int k = 0; k < v.size(); k++)
		{
			transcript &t = v[k];
			if(t2i.find(t.transcript_id) == t2i.end()) continue;
			int x = t2i[t.transcript_id];
			if(items[x].tpm < min_tpm) continue;
			if(items[x].numreads < min_numreads) continue;
			t.TPM = items[x].tpm;
			t.write(cout);
		}
	}
	return 0;
}

int gtfquant::compare()
{
	build_msd();
	vector<double> vx;
	vector<double> vy;

	for(MSDD::iterator it = msd.begin(); it != msd.end(); it++)
	{
		double x = it->second.first + 1.0;
		double y = it->second.second + 1.0;
		x = log(x);
		y = log(y);
		vx.push_back(x);
		vy.push_back(y);
	}

	// compute average
	double ax = 0;
	double ay = 0;
	for(int i = 0; i < vx.size(); i++)
	{
		ax += vx[i];
		ay += vy[i];
	}
	ax /= vx.size();
	ay /= vy.size();

	// compute covariance
	double cov = 0;
	double devx = 0;
	double devy = 0;
	for(int i = 0;i < vx.size(); i++)
	{
		cov += (vx[i] - ax) * (vy[i] - ay);
		devx += (vx[i] - ax) * (vx[i] - ax);
		devy += (vy[i] - ay) * (vy[i] - ay);
	}
	cov /= vx.size();
	devx = sqrt(devx / vx.size());
	devy = sqrt(devy / vy.size());

	double pearson = cov / devx / devy;

	printf("pearson = %.3lf\n", pearson);

	return 0;
}

int gtfquant::build_msd()
{
	msd.clear();
	for(int i = 0; i < gm.genes.size(); i++)
	{
		vector<transcript> &v = gm.genes[i].transcripts;
		for(int k = 0; k < v.size(); k++)
		{
			transcript &t = v[k];
			string s = t.transcript_id;

			if(msd.find(s) != msd.end()) continue;
			double tpm_pred = t.TPM;
			double tpm_quant = 0;

			if(t2i.find(s) != t2i.end())
			{
				int x = t2i[t.transcript_id];
				tpm_quant = items[x].tpm;
			}
			msd.insert(PSDD(s, PDD(tpm_pred, tpm_quant)));
		}
	}

	for(map<string, int>::iterator it = t2i.begin(); it != t2i.end(); it++)
	{
		string s = it->first;
		double tpm_quant = it->second;
		if(msd.find(s) != msd.end()) continue;
		msd.insert(PSDD(s, PDD(0, tpm_quant)));
	}
	return 0;
}
