#include "gtfcuff.h"
#include "genome.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>

using namespace std;

gtfcuff::gtfcuff(const string &cufffile)
{
	read_cuff(cufffile);
	build_cuff_index();
}

int gtfcuff::assign_pred(const string &file)
{
	genome gm(file);
	vpred = gm.collect_transcripts();
	build_pred_index();
	return 0;
}

int gtfcuff::assign_ref(const string &file)
{
	genome gm(file);
	vref = gm.collect_transcripts();
	build_ref_index();
	return 0;
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

int gtfcuff::build_cuff_index()
{
	t2i.clear();
	for(int i = 0; i < items.size(); i++)
	{
		string s = items[i].transcript_id;
		t2i.insert(pair<string, int>(s, i));
	}
	return 0;
}

int gtfcuff::build_pred_index()
{
	t2p.clear();
	for(int i = 0; i < vpred.size(); i++)
	{
		transcript &t = vpred[i];
		string s = t.transcript_id;
		t2p.insert(pair<string, int>(s, i));
	}
	return 0;
}

int gtfcuff::build_ref_index()
{
	t2r.clear();
	for(int i = 0; i < vref.size(); i++)
	{
		transcript &t = vref[i];
		string s = t.transcript_id;
		t2r.insert(pair<string, int>(s, i));
	}
	return 0;
}

int gtfcuff::classify(const string &fn1, const string &fn2)
{
	ofstream f1(fn1.c_str());
	ofstream f2(fn2.c_str());

	for(int k = 0; k < vpred.size(); k++)
	{
		string s = vpred[k].transcript_id;
		bool b = false;
		if(t2i.find(s) != t2i.end())
		{
			if(items[t2i[s]].code == '=') b = true;
			else b = false;
		}
		if(b == true) vpred[k].write(f1);
		else vpred[k].write(f2);
	}

	f1.close();
	f2.close();
	return 0;
}

int gtfcuff::roc(int refsize)
{
	if(items.size() == 0) return 0;

	sort(items.begin(), items.end(), cuffitem_cmp_coverage);

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

int gtfcuff::quant()
{
	vector<double> qpred;
	vector<double> qref;

	int pcnt = 0;
	int rcnt = 0;
	int prcnt = 0;

	set<string> rr;
	for(int i = 0; i < vpred.size(); i++)
	{
		transcript &t = vpred[i];
		double ptpm = t.TPM;
		double rtpm = 0;
		string s = t.transcript_id;
		bool b = false;
		if(t2i.find(s) != t2i.end())
		{
			int j = t2i[s];
			if(items[j].code == '=')
			{
				string r = items[j].ref_transcript_id;
				if(t2r.find(r) != t2r.end())
				{
					int k = t2r[r];
					rtpm = vref[k].TPM;
					rr.insert(r);
					b = true;
					prcnt++;
				}
			}
		}
		if(b == false) pcnt++;
		qpred.push_back(ptpm);
		qref.push_back(rtpm);
	}

	for(int i = 0; i < vref.size(); i++)
	{
		transcript &t = vref[i];
		string s = t.transcript_id;
		if(rr.find(s) != rr.end()) continue;
		double ptpm = 0;
		double rtpm = t.TPM;
		qpred.push_back(ptpm);
		qref.push_back(rtpm);
		rcnt++;
	}
	
	// compute pearson
	vector<double> vx;
	vector<double> vy;

	for(int i = 0; i < qpred.size(); i++)
	{
		double x = log(qpred[i] + 1.0);
		double y = log(qref[i] + 1.0);
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

	printf("pearson = %.3lf, counts of (prediction, common, reference) = (%d, %d, %d)\n", pearson, pcnt, prcnt, rcnt);

	return 0;
}
