/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __VERTEX_INFO__
#define __VERTEX_INFO__

#include <stdint.h>

class vertex_info
{
public:
	vertex_info();
	vertex_info(int l);
	vertex_info(const vertex_info &vi);

public:
	int32_t pos;		// position
	int32_t lpos;		// left position
	int32_t rpos;		// right position
	double stddev;		// standard deviation of read coverage
	int length;			// length of this partial exon
	int sdist;			// shortest distance to s
	int tdist;			// shortest distance to t
	int type;			// for various usage
	char lstrand;		// left side strand
	char rstrand;		// right side strand	
};

#endif
