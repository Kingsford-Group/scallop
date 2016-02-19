#ifndef __INTERVAL_MAP_H__
#define __INTERVAL_MAP_H__

// boost::interval map
#include "boost/icl/interval_map.hpp"
#include "boost/icl/split_interval_map.hpp"
using namespace boost;
using namespace std;

typedef icl::right_open_interval<int32_t> ROI;
typedef icl::split_interval_map<int32_t, int32_t, icl::partial_absorber, less, icl::inplace_plus, icl::inter_section, ROI> interval_map;
typedef interval_map::const_iterator ICI;
typedef pair<ICI, ICI> PICI;

// if p is inside an interval, split this interval into 2
int create_split(interval_map &imap, int32_t p);

// return the overlap at position p
int compute_overlap(const interval_map &imap, int32_t p);

// find the leftmost iterator whose upper posistion <= x
ICI locate_right_iterator(const interval_map &imap, int32_t x);

// find the rightmost interval whose lower position >= x
ICI locate_left_iterator(const interval_map &imap, int32_t x);

// locate boundary iterators
PICI locate_boundary_iterators(const interval_map &imap, int32_t x, int32_t y);

// return the sum of the lengths of intervals from p to q (include q)
int compute_coverage(const interval_map &imap, ICI &p, ICI &q);

// return the maximum overlap of the intervals from p to q (include q)
int compute_max_overlap(const interval_map &imap, ICI &p, ICI &q);

// return the sum of the overlap of the intervals from p to q (include q)
int compute_sum_overlap(const interval_map &imap, ICI &p, ICI &q);


// testing
int test_interval_map();

#endif
