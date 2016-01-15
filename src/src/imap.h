#ifndef __IMAP_H__
#define __IMAP_H__

// boost::interval map
#include <boost/icl/interval_map.hpp>
using namespace boost;
using namespace std;

typedef icl::right_open_interval<int32_t> ROI;
typedef icl::interval_map<int32_t, int32_t, icl::partial_absorber, less, icl::inplace_plus, icl::inter_section, ROI> imap_t;

// return the overlap at position p
int compute_overlap(const imap_t &imap, int32_t p);
// return the sum of overlap in [p,q)
int cumulate_overlap(const imap_t &imap, int32_t p, int32_t q, int32_t t);
// return the max of overlap in [p,q)
int maximum_overlap(const imap_t &imap, int32_t p, int32_t q, int32_t t);
// return the number of positions with positive overlap in [p,q)
int compute_coverage(const imap_t &imap, int32_t p, int32_t q, int32_t t);

#endif
