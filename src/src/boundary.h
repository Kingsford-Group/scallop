#ifndef __POSITION_H__
#define __POSITION_H__

#include <stdint.h>

class boundary
{
public:
	boundary();
	boundary(int _t, int32_t _p);
	boundary(int _t, int32_t _p, int32_t _c, uint32_t _min, uint32_t _max);
	boundary(const boundary &p);

public:
	int type;			// type
	int32_t pos;		// position
	int32_t count;		// number of hits having this splice boundary
	uint32_t min_qual;	// minimum quality
	uint32_t max_qual;	// maximum quality

	int32_t total;		// total number of reads in [lpos - rpos]
	int32_t lpos;		// left position of background interval
	int32_t rpos;		// right position of background interval

public:
	int print();
};

#endif
