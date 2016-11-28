#include "cuffroc.h"
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

cuffroc::cuffroc(const string &file, int n)
{
	refsize = n;
	read(file);
}

int cuffroc::read(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail()) return -1;

	char line[10240];

	while(fin.getline(line, 10240, '\n'))
	{
		item x(line);
		if(x.code == '@') continue;
		items.push_back(x);
	}
	return 0;
}

int cuffroc::solve()
{
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
