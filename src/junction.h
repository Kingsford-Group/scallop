/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __BRIDGE_H__
#define __BRIDGE_H__

#include <stdint.h>
#include <string>

using namespace std;

class junction
{
public:
	junction();
	junction(int64_t _p);
	junction(int64_t _p, int32_t _c);
	junction(const junction &p);

	bool operator<(const junction &x) const;

public:
	int32_t lpos;		// left position [left, right)
	int32_t rpos;		// right position
	int count;			// number of hits having this splice junction
	char strand;		// strandness of this junction
	int nm;				// total mismatch

	int lexon;			// region index corresponds to lpos
	int rexon;			// region index corresponds to rpos

public:
	int print(const string &chrm, int index) const;
};

bool junction_cmp_length(const junction &x, const junction &y);

#endif
