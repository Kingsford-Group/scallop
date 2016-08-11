#ifndef __VERTEX_INFO__
#define __VERTEX_INFO__

#include <cstdint>

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
	double reliability;	// whether the coverage is reliable
	bool infer;			// whether this vertex has been inferred
	int type;			// for various usage
	bool adjust;		// whether the coverage of this region has been adjusted 
};

#endif
