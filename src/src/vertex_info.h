#ifndef __VERTEX_INFO__
#define __VERTEX_INFO__

class vertex_info
{
public:
	vertex_info();
	vertex_info(int l);
	vertex_info(const vertex_info &vi);

public:
	double stddev;		// standard deviation of read coverage
	int length;			// length of this partial exon
	double scalor1;		// not a rectangle
	double scalor2;		// not a rectangle
	int length1;		// length of left part
	int length2;		// length of right part
};

#endif
