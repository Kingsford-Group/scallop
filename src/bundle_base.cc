/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cassert>
#include <cstdio>
#include <cmath>

#include "bundle_base.h"

bundle_base::bundle_base()
{
	tid = -1;
	chrm = "";
	lpos = INT32_MAX;
	rpos = 0;
	strand = '.';
}

bundle_base::~bundle_base()
{}

int bundle_base::add_hit(const hit &ht)
{
	// store new hit
	hits.push_back(ht);

	// calcuate the boundaries on reference
	if(ht.pos < lpos) lpos = ht.pos;
	if(ht.rpos > rpos) rpos = ht.rpos;

	// set tid
	if(tid == -1) tid = ht.tid;
	assert(tid == ht.tid);

	// set strand
	if(hits.size() <= 1) strand = ht.strand;
	assert(strand == ht.strand);

	// DEBUG
	/*
	if(strand != ht.strand)
	{
		printf("strand = %c, ht.strand = %c, ht.xs = %c,\n", strand, ht.strand, ht.xs);
	}
	*/

	// add intervals
	vector<int64_t> vm;
	vector<int64_t> vi;
	vector<int64_t> vd;
	ht.get_mid_intervals(vm, vi, vd);

	//ht.print();
	for(int k = 0; k < vm.size(); k++)
	{
		int32_t s = high32(vm[k]);
		int32_t t = low32(vm[k]);
		//printf(" add interval %d-%d\n", s, t);
		mmap += make_pair(ROI(s, t), 1);
	}

	for(int k = 0; k < vi.size(); k++)
	{
		int32_t s = high32(vi[k]);
		int32_t t = low32(vi[k]);
		imap += make_pair(ROI(s, t), 1);
	}

	for(int k = 0; k < vd.size(); k++)
	{
		int32_t s = high32(vd[k]);
		int32_t t = low32(vd[k]);
		imap += make_pair(ROI(s, t), 1);
	}

	return 0;
}

bool bundle_base::overlap(const hit &ht) const
{
	if(mmap.find(ROI(ht.pos, ht.pos + 1)) != mmap.end()) return true;
	if(mmap.find(ROI(ht.rpos - 1, ht.rpos)) != mmap.end()) return true;
	return false;
}

int bundle_base::clear()
{
	tid = -1;
	chrm = "";
	lpos = INT32_MAX;
	rpos = 0;
	strand = '.';
	hits.clear();
	mmap.clear();
	imap.clear();
	trsts.clear();
	return 0;
}

vector<int> bundle_base::hits_on_transcripts(vector <uint32_t> &pair_mapping){
	int mapped = 0;
  float total_distance = 0;
  float norm_distance = 0;
  int mapped_pairs = 0;
  int pairs = 0;
  //uint64_t max_value = 0;
  //vector<uint64_t> read_index;
	//cout << "\n\n";
  for(int i=0; i<hits.size(); i++){
    //read_index.push_back(stoul(hits[i].qname.substr(hits[i].qname.find(".")+1)));
    //if(read_index[read_index.size()-1] > max_value) max_value = read_index[read_index.size()-1];
    //cout << read_index[read_index.size()-1] << "\t" << max_value << endl;

		for(int j=0; j<trsts.size(); j++){
			if(hits[i].maps_to_transcript(trsts[j])){
        total_distance += hits[i].nm;
        norm_distance += 1.0 * hits[i].nm / hits[i].qlen;
				mapped++;
				break;
      }
		}
	}
  //vector<uint32_t> pair_mapping(max_value,-1);
  for(int i=0; i<hits.size(); i++){
    uint32_t read_index = stoul(hits[i].qname.substr(hits[i].qname.find(".")+1));
    assert(read_index < pair_mapping.size());
    if(pair_mapping[read_index] != -1){
      int ip = pair_mapping[read_index];
      pairs++;
      for(int j=0; j<trsts.size(); j++){
        if(hits[i].maps_to_transcript(trsts[j]) && hits[ip].maps_to_transcript(trsts[j])){
          mapped_pairs++;
          break;          
        }
      }
    }
    pair_mapping[read_index] = i;
  }

    for(int i=0; i<hits.size(); i++){
      uint32_t read_index = stoul(hits[i].qname.substr(hits[i].qname.find(".")+1));
      pair_mapping[read_index] = -1;
    }

  /*for(int i=0; i<hits.size()-1; i++){
    //for(int ip=i+1; ip<=i+1; ip++){
    for(int ip=i+1; ip<hits.size(); ip++){
      if(config::verbose >= 2){
        cerr << i << ": " << hits[i].qname << "=?" << hits[ip].qname << " && " << hits[i].hi << "=?" << hits[ip].hi << endl;
      }
      if(hits[i].qname == hits[ip].qname){
      //if(hits[i].qname == hits[ip].qname && hits[i].hi == hits[ip].hi){
        pairs++;
        if(ip != i+1){
          //cout << "Found one thats not i+1 " << i << "," << ip <<  "!!"  << hits[i].qname <<" == " << hits[ip].qname << " " << hits[i].hi << " == " <<hits[ip].hi<< endl;
        }
        for(int j=0; j<trsts.size(); j++){
          if(hits[i].maps_to_transcript(trsts[j]) && hits[ip].maps_to_transcript(trsts[j])){
            mapped_pairs++;
            break;
          }
        }
        break;
      }
    }
  }*/
	vector<int> rtn;
	rtn.push_back(mapped);
  rtn.push_back(total_distance);
  rtn.push_back(norm_distance);
  rtn.push_back(mapped_pairs);
  rtn.push_back(pairs);
	return rtn;
}
