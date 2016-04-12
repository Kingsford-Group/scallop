#ifndef __BUNDLE_EXTEND_H__
#define __BUNDLE_EXTEND_H__

#include <fstream>

#include "bundle.h"
#include "genome.h"

using namespace std;

class bundle_extend : public bundle
{
public:
	bundle_extend(const bundle_base &b, const genome &g);

public:
	const genome &gn;

public:
	int build();

};

#endif
