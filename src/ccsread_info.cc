#include "ccsread_info.h"

#include <iostream>
#include <sstream>

ccsread_info::ccsread_info()
{}

ccsread_info::ccsread_info(const string &line)
{
	parse(line);
}

int ccsread_info::parse(const string &line)
{
	if(line.size() == 0) return 0;
	if(line[0] != '@') return 0;

	stringstream sstr(line.substr(1, line.size() - 1));

	getline(sstr, qname, ' ');

	string s;
	while(getline(sstr, s, ';'))
	{
		stringstream ss(s);
		string s1, s2;
		getline(ss, s1, '=');
		getline(ss, s2, '=');

		if(s1 == "strand") strand = s2[0];
		else if(s1 == "fiveseen") fiveseen = (s2 == "0") ? false : true;
		else if(s1 == "polyAseen") polyAseen = (s2 == "0") ? false : true;
		else if(s1 == "threeseen") threeseen = (s2 == "0") ? false : true;
		else if(s1 == "fiveend") fiveend = atoi(s2.c_str());
		else if(s1 == "polyAend") polyAend = atoi(s2.c_str());
		else if(s1 == "threeend") threeend = atoi(s2.c_str());
		else if(s1 == "primer") primer = atoi(s2.c_str());
		else if(s1 == "chimera") chimera = (s2 == "0") ? false : true;
	}

	return 0;
}

int ccsread_info::print() const
{
	printf("%s ", qname.c_str());
	printf("strand=%c;", strand);
	printf("fiveseen=%c;", fiveseen ? '1' : '0');
	printf("polyAseen=%c;", polyAseen ? '1' : '0');
	printf("threeseen=%c;", threeseen ? '1' : '0');
	printf("fiveend=%d;", fiveend);
	printf("polyAend=%d;", polyAend);
	printf("threeend=%d;", threeend);
	printf("primer=%d;", primer);
	printf("chimera=%c", chimera ? '1' : '0');
	printf("\n");
	return 0;
}
