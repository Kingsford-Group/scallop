#include "super_hit.h"

int super_hit::add_hit(int x)
{
	hit_list.insert(x);
	return 0;
}
