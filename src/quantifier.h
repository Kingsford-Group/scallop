/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __QUANTIFIER_H__
#define __QUANTIFIER_H__

#include "hit.h"
#include "path.h"

class quantifier
{
public:
	quantifier(const vector<hit> &v, vector<path> &p);

public:
	const vector<hit> &hits;	// reference to hits
	vector<path> &paths;		// paths to be quantified
};

#endif
