#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__

#include <vector>

using namespace std;


#define SMIN 0.000000001

typedef vector< vector<double> > VVD;

class algebra
{
public:
	algebra(const VVD &m);
	~algebra();

public:
	VVD mat;

public:
	int partial_eliminate();
	int full_eliminate();
	int print();

private:
	int partial_eliminate(int r, int l);
	int full_eliminate(int r, int l);
	int choose_main_element(int r, int l);
	int choose_main_row(int r, int l);
	int normalize_row(int r, int l);
	int add_to_row(int s, int t, int l, double c);
	int exchange_row(int s, int t);
	int exchange_column(int s, int t);
};

#endif
