#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "config.h"
#include "scallop.h"

int main(int argc, char **argv)
{
	if(argc != 3) return 0;

	config *conf = new config(argv[1]);

	scallop sc(conf);
	sc.process(argv[2]);

	delete conf;

    return 0;
}
