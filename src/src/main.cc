#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include "manager.h"
#include "config.h"

using namespace std;

int main(int argc, const char **argv)
{
	bool b = parse_arguments(argc, argv);
	if(b == false) return 0;

	manager sc;
	sc.process();

    return 0;
}
