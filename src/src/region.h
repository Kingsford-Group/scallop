#ifndef __REGION_H__
#define __REGION_H__

#include <stdint.h>
#include <vector>
#include <boost/icl/interval_map.hpp>

using namespace std;
using namespace boost;

typedef pair<size_t, size_t> PT;
typedef icl::right_open_interval<int32_t> ROI;
typedef icl::interval_map<int32_t, int32_t, icl::partial_absorber, less, icl::inplace_plus, icl::inter_section, ROI> imap_t;

class region
{
public:
	region(int32_t _lpos, int32_t _rpos, const imap_t *_imap);
	region(const region &r);
	region& operator=(const region &r);
	~region();

public:
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	const imap_t *imap;				// pointer to a interval map
};

#endif
