#include "region.h"

region::region(int32_t _lpos, int32_t _rpos, const imap_t *_imap)
	:lpos(_lpos), rpos(_rpos), imap(_imap)
{}

region::region(const region &r)
	:lpos(r.lpos), rpos(r.rpos), imap(r.imap)
{}

region& region::operator=(const region &r)
{
	lpos = r.lpos;
	rpos = r.rpos;
	imap = r.imap;
	return *this;
}

region::~region()
{}
