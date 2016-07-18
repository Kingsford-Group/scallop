#include "partial_exon.h"
#include "util.h"

partial_exon::partial_exon(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype)
	: lpos(_lpos), rpos(_rpos), ltype(_ltype), rtype(_rtype)
{}

string partial_exon::label() const
{
	string l = tostring((lpos + 1) % 100000);
	string r = tostring(rpos % 100000);
	return (l + "-" + r);
}

int partial_exon::print(int index) const
{
	printf("partial_exon %d: [%d-%d), type = (%d, %d), length = %d, ave-abd = %.1lf, std-abd = %.1lf\n",
			index, lpos, rpos, ltype, rtype, rpos - lpos, ave_abd, dev_abd);
	return 0;
}
