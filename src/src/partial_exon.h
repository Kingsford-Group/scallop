/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __PARTIAL_EXON_H__
#define __PARTIAL_EXON_H__

#include <stdint.h>
#include <vector>
#include <string>

using namespace std;

class partial_exon
{
public:
	partial_exon(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype);

public:
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	int ltype;						// type of the left boundary
	int rtype;						// type of the right boundary

	double ave;						// average abundance
	double dev;						// standard-deviation of abundance

public:
	string label() const;
	int print(int index) const;
};

#endif
