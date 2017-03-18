/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "matrix.h"

vector<double> solve_linear_system(const vector< vector<double> > &A, const vector<double> &B)
{
	int n = A.size();
	ublas::matrix<double> AA(n, n);
	for(int i = 0; i < n; i++)
	{
		assert(A[i].size() == n);
		for(int j = 0; j < n; j++)
		{
			AA(i, j) = A[i][j];
		}
	}

	assert(B.size() == n);
	ublas::vector<double> BB(n);
	for(int i = 0; i < n; i++)
	{
		BB(i) = B[i];
	}

	ublas::permutation_matrix<int> PM(n);
	ublas::lu_factorize(AA, PM);
	ublas::lu_substitute(AA, PM, BB);

	vector<double> X;
	for(int i = 0; i < n; i++)
	{
		X.push_back(BB[i]);
	}
	return X;
}

int test_matrix()
{
	ublas::matrix<double> A(2, 2);
	A(0, 0) = A(1, 1) = 2;
	A(0, 1) = A(1, 0) = 1;
	cout << A << endl;

	ublas::vector<double> b(2);
	b(0) = 20;
	b(1) = 18;
	cout << b << endl;

	ublas::permutation_matrix<int> pm(2);
	ublas::lu_factorize(A, pm);
	ublas::lu_substitute(A, pm, b);

	cout << A << endl;
	cout << b << endl;

	return 0;
}
