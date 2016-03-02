#ifndef __STRINGTIE_H__
#define __STRINGTIE_H__

#include "assembler.h"

class stringtie : public assembler
{
public:
	stringtie(splice_graph &g);
	virtual ~stringtie();

public:
	int assemble();

private:
	int greedy();
};

#endif
