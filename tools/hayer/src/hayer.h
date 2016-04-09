#ifndef __HAYER_H__
#define __HAYER_H__

#include "genome.h"

#include <string>
#include <map>

using namespace std;

typedef map<string, int> MSI;
typedef pair<string, int> PSI;

class hayer
{
public:
	int process(const string &file, const string &s);
	int write(const string &file);

public:
	int read_genome(const string &file);
	int read_transcripts(const string &file);
	int build_genome();
	string parse_transcript_id(const string &s);
	int parse_location(const string &s, string &chrm, int &start, int &end);

private:
	vector<transcript> transcripts;
	genome gm;
};

#endif
