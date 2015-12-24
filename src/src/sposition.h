#ifndef __SPOSITION_H__
#define __SPOSITION_H__

#include <stdint.h>

class sposition
{
public:
	int32_t pos;		// splice position
	int32_t count;		// number of hits having this splice position
	uint32_t min_qual;	// minimum quality
	uint32_t max_qual;	// maximum quality

public:
	sposition();
	sposition(int32_t _p, int32_t _c, uint32_t _min, uint32_t _max);
	sposition(const sposition &sp);

public:	
	int print();
};

#endif
