#include "bundle_extend.h"

bundle_extend::bundle_extend(const bundle_base &b, const genome &g)
	: gn(g), bundle(b)
{}

int bundle_extend::build()
{
	return 0;
}
