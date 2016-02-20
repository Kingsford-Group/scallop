#include "gtf_gene.h"

int gtf_gene::add_exon(const gtf_exon &ge)
{
	exons.push_back(ge);
	return 0;
}
