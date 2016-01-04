#ifndef __BRIDGE_H__
#define __BRIDGE_H__

#include <stdint.h>

class bridge
{
public:
	bridge();
	bridge(int _t, int64_t _p);
	bridge(int _t, int64_t _p, int32_t _c, uint32_t _min, uint32_t _max);
	bridge(const bridge &p);

	bool operator<(const bridge &x) const;

public:
	int type;			// type
	int32_t lpos;		// left position [left, right)
	int32_t rpos;		// right position
	int32_t count;		// number of hits having this splice bridge
	uint32_t min_qual;	// minimum quality
	uint32_t max_qual;	// maximum quality
	uint32_t score;		// score

public:
	int print();
};

#endif
