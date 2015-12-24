#ifndef __POSITION_H__
#define __POSITION_H__

#include <stdint.h>

class position
{
public:
	position(int32_t _p);
	position(const position &p);
	virtual ~position();

public:
	int32_t pos;

public:
	virtual int print() = 0;
};

class splice_pos: public position
{
public:
	int32_t count;		// number of hits having this splice position
	uint32_t min_qual;	// minimum quality
	uint32_t max_qual;	// maximum quality

public:
	splice_pos(int32_t _p);
	splice_pos(int32_t _p, int32_t _c, uint32_t _min, uint32_t _max);
	splice_pos(const splice_pos &sp);

public:	
	int print();
};

#endif
