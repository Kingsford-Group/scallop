#ifndef __CONFIG_H__
#define __CONFIG_H__


#include <stdint.h>
using namespace std;

class config
{
public:
	config(const char * conf_file);
	int load_config(const char * conf_file);

public:
	int32_t min_bundle_gap;
};

#endif
