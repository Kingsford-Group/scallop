#ifndef __SHANNON_H__
#define __SHANNON_H__

#include "assembler.h"

class shannon : public assembler
{
public:
	shannon(splice_graph &g);
	virtual ~shannon();

public:
	int assemble();
};

#endif
