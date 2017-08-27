#include "metaassembler.h"
#include "previewer.h"
#include "sgraph_compare.h"

int metaassembler::preassemble()
{
	previewer pv1(input_file1);
	pv1.preview();
	library_type1 = library_type;

	previewer pv2(input_file2);
	pv2.preview();
	library_type2 = library_type;

	library_type = library_type1;
	input_file = input_file1;
	assembler asmb1;
	asmb1.preassemble();
	grlist1 = asmb1.grlist;
	hslist1 = asmb1.hslist;

	library_type = library_type2;
	input_file = input_file2;
	assembler asmb2;
	asmb2.preassemble();
	grlist2 = asmb2.grlist;
	hslist2 = asmb2.hslist;

	print();

	return 0;
}

int metaassembler::assemble()
{
	compare_boundaries();
	return 0;
}

int metaassembler::postassemble()
{
	return 0;
}

int metaassembler::compare_boundaries()
{
	vector< set<int32_t> > brlist1;
	vector< set<int32_t> > brlist2;

	for(int k = 0; k < grlist1.size(); k++)
	{
		set<int32_t> s = grlist1[k].get_boundaries();
		brlist1.push_back(s);
	}

	for(int k = 0; k < grlist2.size(); k++)
	{
		set<int32_t> s = grlist2[k].get_boundaries();
		brlist2.push_back(s);
	}

	for(int i = 0; i < brlist1.size(); i++)
	{
		int index = -1;
		int maxbr = -1;
		set<int32_t> &s1 = brlist1[i];
		vector<int32_t> v(s1.size());
		for(int j = 0; j < brlist2.size(); j++)
		{
			set<int32_t> &s2 = brlist2[j];	
			vector<int32_t>::iterator it = set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), v.begin());
			int m = it - v.begin();

			if(m > maxbr)
			{
				maxbr = m;
				index = j;
			}
		}

		printf("intersection between graph %d (%lu boundaries) and graph %d (%lu boundaries) = %d boundaries\n", 
				i, s1.size(), index, brlist2[index].size(), maxbr);
	}

	sgraph_compare sc(grlist1[593], grlist2[638]);
	sc.compare("meta.593.639.tex");

	return 0;
}

int metaassembler::print()
{
	printf("bamfiles contain %lu and %lu splice graphs\n", grlist1.size(), grlist2.size());
	return 0;
}
