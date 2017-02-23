#include "gtfcuff.h"
#include "genome.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <cfloat>

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

int gtfcuff::read_quant(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail()) return -1;

	char line[10240];
	while(fin.getline(line, 10240, '\n'))
	{
		quantitem x(line);
		qitems.push_back(x);
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

int gtfcuff::build_quant_index()
{
	t2q.clear();
	for(int i = 0; i < qitems.size(); i++)
	{
		string s = qitems[i].transcript_id;
		//printf("%s -> %d\n", s.c_str(), i);
		t2q.insert(pair<string, int>(s, i));
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

int gtfcuff::split(const string &fn1, const string &fn2)
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
	return roc_trunc(refsize, -1, DBL_MAX);
}

int gtfcuff::roc_trunc(int refsize, double min_coverage, double max_coverage)
{
	if(items.size() == 0) return 0;

	vector<cuffitem> vt;
	for(int i = 0; i < items.size(); i++)
	{
		if(items[i].coverage < min_coverage) continue;
		if(items[i].coverage > max_coverage) continue;
		vt.push_back(items[i]);
	}

	sort(vt.begin(), vt.end(), cuffitem_cmp_coverage);

	int correct = 0;
	for(int i = 0; i < vt.size(); i++) if(vt[i].code == '=') correct++;

	double max_sen = 0;
	double max_pre = 0;
	double max_cov = 0;
	int max_len = 0;
	int max_correct = 0;
	int max_size = 0;
	double sen0 = correct * 100.0 / refsize;
	for(int i = 0; i < vt.size(); i++)
	{
		double sen = correct * 100.0 / refsize;
		double pre = correct * 100.0 / (vt.size() - i);

		if(sen * 10.0 < sen0) break;

		if(i % 100 == 0)
		{
			printf("ROC: reference = %d prediction = %lu correct = %d sensitivity = %.2lf precision = %.2lf | coverage = %.3lf, length = %d\n",
				refsize, vt.size() - i, correct, sen, pre, vt[i].coverage, vt[i].length);
		}

		if(vt[i].code == '=') correct--;
	}

	return 0;
}

int gtfcuff::roc_quant(const string &qfile, double min_tpm, double max_tpm)
{
	read_quant(qfile);
	vector<quantitem> vq;
	for(int i = 0; i < qitems.size(); i++)
	{
		if(qitems[i].tpm < min_tpm) continue;
		if(qitems[i].tpm > max_tpm) continue;
		vq.push_back(qitems[i]);
	}
	qitems = vq;
	build_quant_index();

	if(items.size() == 0) return 0;
	if(qitems.size() == 0) return 0;

	int refsize = qitems.size();
	vector<cuffitem> vt;
	for(int i = 0; i < items.size(); i++)
	{
		cuffitem ci = items[i];
		string s = ci.ref_transcript_id;
		if(t2q.find(s) == t2q.end()) ci.code = 'x';
		vt.push_back(ci);
	}

	sort(vt.begin(), vt.end(), cuffitem_cmp_coverage);

	int correct = 0;
	for(int i = 0; i < vt.size(); i++) if(vt[i].code == '=') correct++;

	double max_sen = 0;
	double max_pre = 0;
	double max_cov = 0;
	int max_len = 0;
	int max_correct = 0;
	int max_size = 0;
	double sen0 = correct * 100.0 / refsize;
	for(int i = 0; i < vt.size(); i++)
	{
		double sen = correct * 100.0 / refsize;
		double pre = correct * 100.0 / (vt.size() - i);

		if(sen * 10.0 < sen0) break;

		if(i % 100 == 0)
		{
			printf("ROC: reference = %d prediction = %lu correct = %d sensitivity = %.2lf precision = %.2lf | coverage = %.3lf, length = %d\n",
				refsize, vt.size() - i, correct, sen, pre, vt[i].coverage, vt[i].length);
		}

		if(vt[i].code == '=') correct--;
	}

	return 0;
}

int gtfcuff::acc(int refsize)
{
	if(items.size() == 0) return 0;
	sort(items.begin(), items.end(), cuffitem_cmp_coverage);
	int n = items.size() / 3;

	int correct1 = 0;
	int correct2 = 0;
	int correct3 = 0;
	for(int i = 0; i < items.size(); i++)
	{
		if(i >= n * 0 && i < n * 1 && items[i].code == '=') correct1++;
		if(i >= n * 1 && i < n * 2 && items[i].code == '=') correct2++;
		if(i >= n * 2 && i < n * 3 && items[i].code == '=') correct3++;
	}

	double sen1 = 100.0 * correct1 / refsize;
	double sen2 = 100.0 * correct2 / refsize;
	double sen3 = 100.0 * correct3 / refsize;
	double pre1 = 100.0 * correct1 / n;
	double pre2 = 100.0 * correct2 / n;
	double pre3 = 100.0 * correct3 / n;
	printf("ROC1: reference = %d prediction = %d correct = %d sensitivity = %.2lf precision = %.2lf\n", refsize, n, correct1, sen1, pre1);
	printf("ROC2: reference = %d prediction = %d correct = %d sensitivity = %.2lf precision = %.2lf\n", refsize, n, correct2, sen2, pre2);
	printf("ROC3: reference = %d prediction = %d correct = %d sensitivity = %.2lf precision = %.2lf\n", refsize, n, correct3, sen3, pre3);
	return 0;
}

int gtfcuff::acc_quant(const string &qfile, double tpm_threshold)
{
	read_quant(qfile);
	vector<quantitem> vq;
	for(int i = 0; i < qitems.size(); i++)
	{
		if(qitems[i].tpm < tpm_threshold) continue;
		vq.push_back(qitems[i]);
	}

	sort(vq.begin(), vq.end());

	int n = vq.size() / 3;
	set<string> s1;
	set<string> s2;
	set<string> s3;
	set<string> s4;
	for(int i = n * 0; i < n * 1; i++) s1.insert(vq[i].transcript_id);
	for(int i = n * 1; i < n * 2; i++) s2.insert(vq[i].transcript_id);
	for(int i = n * 2; i < n * 3; i++) s3.insert(vq[i].transcript_id);
	//for(int i = n * 3; i < n * 4; i++) s4.insert(vq[i].transcript_id);

	if(items.size() == 0) return 0;

	int refsize = n;
	int correct1 = 0;
	int correct2 = 0;
	int correct3 = 0;
	int correct4 = 0;
	for(int i = 0; i < items.size(); i++)
	{
		cuffitem c = items[i];
		string s = c.ref_transcript_id;
		if(s1.find(s) != s1.end() && c.code == '=') correct1++;
		if(s2.find(s) != s2.end() && c.code == '=') correct2++;
		if(s3.find(s) != s3.end() && c.code == '=') correct3++;
		if(s4.find(s) != s4.end() && c.code == '=') correct4++;
	}

	double sen1 = 100.0 * correct1 / items.size();
	double sen2 = 100.0 * correct2 / items.size();
	double sen3 = 100.0 * correct3 / items.size();
	double sen4 = 100.0 * correct4 / items.size();
	double pre1 = 100.0 * correct1 / refsize;
	double pre2 = 100.0 * correct2 / refsize;
	double pre3 = 100.0 * correct3 / refsize;
	double pre4 = 100.0 * correct4 / refsize;
	printf("ROC1: reference = %d prediction = %lu correct = %d sensitivity = %.2lf precision = %.2lf\n", refsize, items.size(), correct1, sen1, pre1);
	printf("ROC2: reference = %d prediction = %lu correct = %d sensitivity = %.2lf precision = %.2lf\n", refsize, items.size(), correct2, sen2, pre2);
	printf("ROC3: reference = %d prediction = %lu correct = %d sensitivity = %.2lf precision = %.2lf\n", refsize, items.size(), correct3, sen3, pre3);
	//printf("ROC4: reference = %d prediction = %lu correct = %d sensitivity = %.2lf precision = %.2lf\n", refsize, items.size(), correct4, sen4, pre4);
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

					printf("pred = %s, ref = %s, pred-quant = %.3lf, ref-quant = %.3lf\n",
							s.c_str(), r.c_str(), ptpm, rtpm);
				}
			}
		}
		if(b == false)
		{
			printf("pred = %s, ref = %s, pred-quant = %.3lf, ref-quant = %.3lf\n", s.c_str(), "N/A", ptpm, rtpm);
			pcnt++;
		}
			
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
		printf("pred = %s, ref = %s, pred-quant = %.3lf, ref-quant = %.3lf\n", "N/A", s.c_str(), ptpm, rtpm);
		rcnt++;
	}
	
	// compute pearson
	double sumx = 0;
	double sumy = 0;
	for(int i = 0; i < qpred.size(); i++)
	{
		sumx += qpred[i];
		sumy += qref[i];
	}

	vector<double> vx;
	vector<double> vy;

	for(int i = 0; i < qpred.size(); i++)
	{
		double x = log(qpred[i] * 1e6 / sumx + 1.0);
		double y = log(qref[i] * 1e6 / sumy + 1.0);
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

int gtfcuff::classify()
{
	int n = 20;
	vector<int> v1(n, 0);		// correct
	vector<int> v2(n, 0);		// total
	for(int k = 0; k < vpred.size(); k++)
	{
		string s = vpred[k].transcript_id;
		bool b = false;
		if(t2i.find(s) != t2i.end())
		{
			if(items[t2i[s]].code == '=') b = true;
			else b = false;
		}
		int e = vpred[k].exons.size();
		assert(e >= 1);
		if(e >= n) e = n;
		v2[e - 1]++;
		if(b == true) v1[e - 1]++;
	}

	for(int k = 0; k < n; k++)
	{
		printf("exons = %d correct = %d total = %d precision = %.2lf\n", k + 1, v1[k], v2[k], v1[k] * 100.0 / v2[k]);
	}
	return 0;
}

