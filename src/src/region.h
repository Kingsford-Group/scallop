#ifndef __REGION_H__
#define __REGION_H__

#include <stdint.h>
#include <vector>
#include "interval_map.h"
#include "slope.h"

using namespace std;

class region
{
public:
	region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const split_interval_map *_imap);
	~region();

public:
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	int ltype;						// type of the left boundary
	int rtype;						// type of the right boundary
	const split_interval_map *imap;	// pointer to a interval map

	int32_t lcore;					// left core position
	int32_t rcore;					// right core position
	bool empty;						// whether this region is completely spliced

	double ave_abd;					// average abundance
	double dev_abd;					// standard-deviation of abundance

private:
	vector<int> bins;			// average abundance for bins
	vector<slope> slopes;		// slopes

public:
	int build();
	bool left_break() const;
	bool right_break() const;
	string label() const;
	int print(int index) const;
	int print_boundaries(int index) const;

private:
	int init_core_empty();
	int estimate_abundance();
	int estimate_abundance(int ll, int rr, double &ave, double &dev);
	double compute_deviation(const split_interval_map &sim);

	// identify boundaries inside this region
	int compute_bin_abundances();
	int compute_slopes();
	int select_slopes();
};

#endif
