#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "scallop.h"
#include "config.h"

#include <boost/icl/interval_map.hpp>

using namespace boost::icl;
using namespace std;

/*
int test_interval_map()
{
	typedef interval_map<int, int, partial_absorber, less, inplace_plus, inter_section, closed_interval<int> > imap_t;
	imap_t imap;
	imap += make_pair(closed_interval<int>(1, 2), 1);
	imap += make_pair(closed_interval<int>(2, 4), 1);
	imap_t::const_iterator it;
	it = imap.find(1); printf("[%d,%d] => %d\n", it->first.lower(), it->first.upper(), it->second);
	it = imap.find(2); printf("[%d,%d] => %d\n", it->first.lower(), it->first.upper(), it->second);
	it = imap.find(3); printf("[%d,%d] => %d\n", it->first.lower(), it->first.upper(), it->second);
	it = imap.find(4); printf("[%d,%d] => %d\n", it->first.lower(), it->first.upper(), it->second);
	return 0;
}
*/

int main(int argc, char **argv)
{
	//test_interval_map();
	return 0;

	if(argc != 3) return 0;

	load_config(argv[1]);

	scallop sc;
	sc.process(argv[2]);

    return 0;
}
