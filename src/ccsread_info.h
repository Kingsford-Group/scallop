#ifndef CCSREAD_INFO
#define CCSREAD_INFO

#include <string>

using namespace std;

class ccsread_info
{
public:
	ccsread_info(const string &line);

public:
	string qname;
	char strand;
	bool fiveseen;
	bool polyAseen;
	bool threeseen;
	int fiveend;
	int polyAend;
	int threeend;
	int primer;
	bool chimera;

public:
	int parse(const string &line);
	int print() const;
};

#endif
