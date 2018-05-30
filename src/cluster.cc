/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "cluster.h"
#include "config.h"
#include <cassert>
#include <algorithm>

cluster::cluster(const vector<transcript> &v)
	:trs(v)
{}

