#ifndef __COMMON_H__
#define __COMMON_H__

// macros
#define pack(x, y) (int64_t)((((int64_t)(x)) << 32) | ((int64_t)(y)))
#define high32(x) (int32_t)((x) >> 32)
#define low32(x) (int32_t)(((x) << 32) >> 32)


#include <boost/icl/interval_map.hpp>
#include <stdint.h>

using namespace std;
using namespace boost;

typedef pair<size_t, size_t> PT;
typedef icl::right_open_interval<int32_t> ROI;
typedef icl::interval_map<int32_t, int32_t, icl::partial_absorber, less, icl::inplace_plus, icl::inter_section, ROI> imap_t;

int compute_overlap(const imap_t &imap, int32_t p);
int cumulate_overlap(const imap_t &imap, int32_t p, int32_t q, int32_t t);

uint32_t compute_binomial_score(int n, double pr, int x);

#endif
