#include "edge_info.h"

edge_info::edge_info()
	: stddev(1.0), length(0)
{}

edge_info::edge_info(int l)
	: length(l)
{}

edge_info::edge_info(const edge_info &ei)
{
	stddev = ei.stddev;
	length = ei.length;
}
