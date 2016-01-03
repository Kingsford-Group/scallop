#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "region.h"
#include "scallop.h"
#include "config.h"

using namespace std;

int test_interval_map()
{
	imap_t imap;
	imap += make_pair(ROI(1, 4), 2);
	imap += make_pair(ROI(2, 5), 2);
	imap -= make_pair(ROI(2, 3), 1);

	int p = 6;
	imap_t::const_iterator it;
	for(int i = 0; i < p; i++)
	{
		it = imap.find(i);
		if(it != imap.end())
		{
			printf("overlap at %d => %d\n", i, it->second);
		}
	}

	for(it = imap.begin(); it != imap.end(); it++)
	{
		printf("(%d,%d) => %d\n", it->first.lower(), it->first.upper(), it->second);
	}
	return 0;
}

int main(int argc, char **argv)
{
	if(argc != 3) return 0;

	load_config(argv[1]);

	scallop sc;
	sc.process(argv[2]);

    return 0;
}
