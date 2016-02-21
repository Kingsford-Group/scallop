#include "shannon.h"

shannon::shannon(splice_graph &gr)
	: assembler(gr)
{}

shannon::~shannon()
{}

int shannon::assemble()
{
	update_weights();
	return 0;
}

