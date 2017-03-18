/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __BRIDGE_H__
#define __BRIDGE_H__

#include <stdint.h>

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

	int lexon;			// region index corresponds to lpos
	int rexon;			// region index corresponds to rpos

public:
	int print(int index) const;
};

#endif
