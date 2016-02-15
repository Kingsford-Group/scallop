#ifndef __BOUNDARY_H__
#define __BOUNDARY_H__

#include <stdint.h>

class boundary
{
public:
	boundary();
	boundary(int _t, int32_t _p);
	boundary(int _t, int32_t _p, int32_t _c, uint32_t _min, uint32_t _max);
	boundary(const boundary &p);

	bool operator<(const boundary &x) const;

public:
	int type;			// type
	int32_t pos;		// position
	int32_t count;		// number of hits having this splice boundary
	uint32_t min_qual;	// minimum quality
	uint32_t max_qual;	// maximum quality
	uint32_t score;		// might be transfomred from pvalue (for start,end-boundary)

public:
	int print(int index) const;
};

#endif
