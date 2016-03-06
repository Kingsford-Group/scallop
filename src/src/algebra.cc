#include "algebra.h"

#include <cstdio>
#include <cfloat>
#include <cmath>
#include <cassert>

algebra::algebra(const VVD &m)
	:mat(m)
{}

algebra::~algebra()
{}

int algebra::full_eliminate()
{
	for(int i = 0; i < mat.size(); i++)
	{
		int f = full_eliminate(i, i);
		if(f == -1) return i;
	}
	return mat.size();
}

int algebra::partial_eliminate()
{
	int i = 0;
	int j = 0;
	int rank = 0;
	while(i < mat.size() && j < mat[0].size())
	{
		int f = partial_eliminate(i, j);
		if(f == -1) j++;
		else
		{
			rank++;
			i++;
			j++;
		}
	}
	return rank;
}

int algebra::full_eliminate(int r, int l)
{
	int f = choose_main_element(r, l);
	if(f == -1) return -1;
	normalize_row(r, l);
	for(int k = 0; k < mat.size(); k++)
	{
		if(k == r) continue;
		if(fabs(mat[k][l]) < SMIN) continue;
		add_to_row(r, k, l, 0 - mat[k][l] / mat[r][l]);
	}
	return 0;
}

int algebra::partial_eliminate(int r, int l)
{
	int f = choose_main_row(r, l);
	if(f == -1) return -1;
	for(int k = r + 1; k < mat.size(); k++)
	{
		if(fabs(mat[k][l]) < SMIN) continue;
		add_to_row(r, k, l, 0 - mat[k][l] / mat[r][l]);
	}
	return 0;
}

int algebra::add_to_row(int s, int t, int l, double c)
{
	for(int i = l; i < mat[0].size(); i++)
	{
		mat[t][i] += mat[s][i] * c;
	}
	return 0;
}

int algebra::normalize_row(int r, int l)
{
	assert(fabs(mat[r][l]) >= SMIN);
	for(int i = 0; i < l; i++) assert(fabs(mat[r][i]) < SMIN);

	for(int i = l + 1; i < mat[0].size(); i++)
	{
		mat[r][i] = mat[r][i] / mat[r][l];
	}
	mat[r][l] = 1;
	return 0;
}

int algebra::choose_main_row(int r, int l)
{
	if(fabs(mat[r][l]) > SMIN) return 0;
	int k = -1;
	for(int i = r + 1; i < mat.size(); i++)
	{
		if(fabs(mat[i][l]) > SMIN)
		{
			k = i;
			break;
		}
	}
	if(k == -1) return -1;
	exchange_row(r, k);
	return 0;
}

int algebra::choose_main_element(int r, int l)
{
	int ii = -1;
	int jj = -1;
	for(int j = l; j < mat[0].size(); j++)
	{
		for(int i = r; i < mat.size(); i++)
		{
			if(fabs(mat[i][j]) > SMIN)
			{
				ii = i;
				jj = j;
				break;
			}
		}
	}
	if(ii == -1) return -1;
	if(ii != r) exchange_row(r, ii);
	if(jj != l) exchange_column(l, jj);
	return 0;
}

int algebra::exchange_row(int s, int t)
{
	vector<double> v = mat[s];
	mat[s] = mat[t];
	mat[t] = v;
	return 0;
}

int algebra::exchange_column(int s, int t)
{
	vector<double> v;
	for(int i = 0; i < mat.size(); i++) v.push_back(mat[i][s]);
	for(int i = 0; i < mat.size(); i++) mat[i][s] = mat[i][t];
	for(int i = 0; i < mat.size(); i++) mat[i][t] = v[i];
	return 0;
}

int algebra::print_matrix(const VVD &m)
{
	printf("     ");
	for(int j = 0; j < m[0].size(); j++)
	{
		printf("[%2d] ", j + 1);
	}
	printf("\n");

	for(int i = 0; i < m.size(); i++)
	{
		printf("[%2d] ", i + 1);
		for(int j = 0; j < m[i].size(); j++)
		{
			if(fabs(m[i][j]) < SMIN) printf("%4.1lf ", 0.0);
			else printf("%4.1lf ", m[i][j]);
		}
		printf("\n");
	}
	return 0;
}
