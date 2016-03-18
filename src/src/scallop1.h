#ifndef __SCALLOP1_H__
#define __SCALLOP1_H__

#include "assembler.h"

// algorithm: using backward edges
class scallop1 : public assembler
{
public:
	scallop1(splice_graph &gr);
	virtual ~scallop1();

public:
	int assemble();

private:
	int iterate();
	int resolve(const path &px, const path &py, path &qx, path &qy);
};

#endif
