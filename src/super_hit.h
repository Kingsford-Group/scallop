/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __SUPER_HIT_H__
#define __SUPER_HIT_H__

#include "hit.h"

class super_hit
{
public:
	set<int> hit_list;			// list of hits that have the same phasing path

public:
	int add_hit(int x);
};

#endif
