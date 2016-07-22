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
	//int sdist;			// dist to s (along the heaviest path)
	//int tdist;			// dist to t (along the heaviest path)
};

#endif
