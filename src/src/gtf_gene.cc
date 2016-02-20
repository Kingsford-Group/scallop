#include "gtf_gene.h"

int gtf_gene::add_exon(const gtf_exon &ge)
{
	exons.push_back(ge);
	return 0;
}

int gtf_gene::print()
{
	for(int i = 0; i < exons.size(); i++)
	{
		exons[i].print();
	}
	return 0;
}
