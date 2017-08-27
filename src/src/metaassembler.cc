#include "metaassembler.h"
#include "previewer.h"

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
	return 0;
}

int metaassembler::postassemble()
{
	return 0;
}

int metaassembler::print()
{
	printf("bamfiles contain %lu and %lu splice graphs\n", grlist1.size(), grlist2.size());
	return 0;
}
