#ifndef __VERTEX_INFO__
#define __VERTEX_INFO__

class vertex_info
{
public:
	vertex_info();
	vertex_info(int l);
	vertex_info(const vertex_info &vi);

public:
	double stddev;
	int length;
};

#endif
