/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __CLUSTER_H__
#define __CLUSTER_H__

#include "gene.h"

class cluster
{
public:
	cluster(const vector<transcript> &v);

public:
	vector<transcript> trs;

};

#endif
