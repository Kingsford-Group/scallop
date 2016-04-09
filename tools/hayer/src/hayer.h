#ifndef __HAYER_H__
#define __HAYER_H__

#include "genome.h"

#include <string>

using namespace std;

class hayer
{
public:
	hayer(const string &file);

public:
	int read(const string &file);
	string parse_transcript_id(const string &s);
	int parse_location(const string &s, string &chrm, int &start, int &end);
	int write(const string &file);

public:
	genome gm;
};

#endif
