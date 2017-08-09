/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __PREVIEWER_H__
#define __PREVIEWER_H__

#include "hit.h"

#include <fstream>
#include <string>

using namespace std;

class previewer
{
public:
	previewer();
	~previewer();

private:
	samFile *sfn;
	bam_hdr_t *hdr;
	bam1_t *b1t;

public:
	int preview();
};

#endif
