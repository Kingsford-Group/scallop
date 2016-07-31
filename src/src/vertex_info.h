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
	int sdist;			// shortest distance to s
	int tdist;			// shortest distance to t
	bool infer;			// whether this vertex has been inferred
	double scalor;		// for inferred vertices
	int type;			// for various usage
};

#endif
