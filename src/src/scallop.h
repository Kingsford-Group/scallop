#ifndef __SCALLOP_H__
#define __SCALLOP_H__

#include "assembler.h"

class scallop : public assembler
{
public:
	scallop(splice_graph &gr);
	virtual ~scallop();

public:
	int assemble();

private:
	int iterate();
	int resolve(const path &px, const path &py, path &qx, path &qy) const;
};

#endif
