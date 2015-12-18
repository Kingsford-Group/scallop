#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "test_sam.h"

int main(int argc, char **argv)
{
	if(argc != 2) return 0;

	test(argv[1]);
    return 0;
}
