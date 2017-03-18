/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

using namespace std;
namespace ublas = boost::numeric::ublas;

vector<double> solve_linear_system(const vector< vector<double> > &A, const vector<double> &B);
int test_matrix();

#endif
