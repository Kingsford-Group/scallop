#include "edge_info.h"

edge_info::edge_info()
	: stddev(1.0), length(0)
{
	infer = false;
	type = 0;
	jid = -1;
}

edge_info::edge_info(int l)
	: length(l)
{
	infer = false;
	type = 0;
	jid = -1;
}

edge_info::edge_info(const edge_info &ei)
{
	stddev = ei.stddev;
	length = ei.length;
	infer = ei.infer;
	type = ei.type;
	jid = ei.jid;
}
