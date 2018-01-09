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

vector<uint64_t> bundle_base::hits_on_transcripts(unordered_map <string,int> &pair_mapping){
	int mapped = 0;
  float total_distance = 0;
  float norm_distance = 0;
  int mapped_pairs = 0;
  int pairs = 0;
  int reads = 0;
  //uint64_t max_value = 0;
  //vector<uint64_t> read_index;
	//cout << "\n\n";
  //cout << "In Hits " << pair_mapping["BundleIndex"] << endl;
  unordered_map<string,int> hit_index;
  vector<uint64_t> read_count(rpos - lpos + 1, 0);
  //for(int i=0; i<(rpos - lpos + 1); i++){
  //   read_count[i] = 0;
  //}
  for(int i=0; i<hits.size(); i++){
    //cout << "NM: " << hits[i].nm << endl;
    //read_index.push_back(stoul(hits[i].qname.substr(hits[i].qname.find(".")+1)));
    //if(read_index[read_index.size()-1] > max_value) max_value = read_index[read_index.size()-1];
    //cout << read_index[read_index.size()-1] << "\t" << max_value << endl;

    string name = hits[i].qname;
    name += ((hits[i].flag & 0x40)?".1":".2");
    assert(!((hits[i].flag & 0x40) && (hits[i].flag & 0x80)));
    assert((hits[i].flag & 0x40) || (hits[i].flag & 0x80));
	
    bool found = false;  
    //for(int j=0; j<trsts.size(); j++) hits[i].transcripts.push_back(trsts[j]);
    for(int j=0; j<trsts.size(); j++){
			if(hits[i].maps_to_transcript(trsts[j])){
        cerr << name << "\t" << pair_mapping["BundleIndex"] << "\t" << trsts[j].seqname << "\t" << trsts[j].gene_id << "\t" << trsts[j].transcript_id << endl;
        found = true;
        if(pair_mapping.find(name) == pair_mapping.end()){
          //cerr << name << "\t" << trsts[j].seqname << "\t" << trsts[j].gene_id << "\t" << trsts[j].transcript_id<< endl;
          mapped++;
          total_distance += hits[i].nm;
          norm_distance += 1.0 * hits[i].nm / hits[i].qlen;

          /*********/
          //cout << "Start read counts" << endl;
          for(int l=0; l<=hits[i].spos.size(); l++){
            uint64_t start = (l==0)?hits[i].pos:low32(hits[i].spos[l-1]);
            uint64_t end = (l==hits[i].spos.size())?hits[i].rpos:high32(hits[i].spos[l]);
            for(uint64_t k=start; k<end; k++){
              //cout << "l:" << (rpos - lpos + 1) << "\tk:" << k << "\tlpos:" << lpos << endl;
              read_count[k - lpos]++;
            }
          }
          //cout << "Finished" << endl;
          /**********/

        }
        pair_mapping[name] = 1;
				//mapped++;
				//break;
      }
		}
    string seen_name = "seen";
    seen_name += name;
    if(pair_mapping.find(seen_name) == pair_mapping.end()){
      reads++;
      if(!found) cerr << name << "\t" << pair_mapping["BundleIndex"] << "\t-\t-\t-\n";
    }
    pair_mapping[seen_name] = 1;
  //}
  //vector<uint32_t> pair_mapping(max_value,-1);
  //for(int i=0; i<hits.size(); i++){
    if(hit_index.find(hits[i].qname) != hit_index.end()){
      seen_name = "seen";
      seen_name += hits[i].qname;
      if(pair_mapping.find(seen_name) == pair_mapping.end())
        pairs++;
      pair_mapping[seen_name] = 1;
      int ip = hit_index[hits[i].qname];
      bool go = true;
      for(int j=0; j<trsts.size() && go; j++){ 
        if(hits[i].maps_to_transcript(trsts[j]) && hits[ip].maps_to_transcript(trsts[j])){
          if(pair_mapping.find(hits[i].qname) == pair_mapping.end()){
            mapped_pairs++;
          }
          pair_mapping[hits[i].qname] = 1;
          go = false;
          break;          
        }
      }
    }
    hit_index[hits[i].qname] = i;
    /*
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
    */
  }

    for(int i=0; i<hits.size(); i++){
      //uint32_t read_index = stoul(hits[i].qname.substr(hits[i].qname.find(".")+1));
      //pair_mapping[read_index] = -1;
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

  //cout << "Between Hits Loop and Transcripts Loop" << endl;
  /*********
  cout << "Pre Declare" << (rpos - lpos + 1) << endl;
  vector<uint64_t> transcript_count(rpos - lpos + 1, 0);
  cout << "Start init" << endl; 
  //for(uint64_t i=0; i<=(rpos-lpos); i++){
  //  transcript_count[i] = 0;
  //}
  cout << "Start loop 1" << endl;
  for(int j=0; j<trsts.size(); j++){
    for(int i=0; i<trsts[j].exons.size(); i++){
      for(uint64_t k=trsts[j].exons[i].first; k<=trsts[j].exons[i].second; k++){
        transcript_count[k - lpos] += pow(2,j);
        cout << "j:" << j << endl;
      }
    }
  }
  cout << "Finished loop 1" << endl;
  vector<uint64_t> total_read_count_per_bin((uint64_t)pow(2,trsts.size()+1),0);
  vector<uint64_t> total_column_count_per_bin((uint64_t)pow(2,trsts.size()+1),0);
  //for(uint64_t i=0; i<pow(2,trsts.size()+1); i++){
  //  total_read_count_per_bin[i] = 0;
  //  total_column_count_per_bin[i] = 0;
  //}

  cout << "Finished loop 2" << endl;

  for(uint64_t i=0; i<=(rpos-lpos); i++){
    cout << "transcript_count[i]:" << transcript_count[i] << endl;
    total_read_count_per_bin[transcript_count[i]] += read_count[i];
    total_column_count_per_bin[transcript_count[i]] ++;
  }

  cout << "Finished loop 3" << endl;
  //calculate average
  //calculate variance
  //return weighted variance
  //return number of transcripts
  *****/

  vector<bool> in_transcript(trsts.size(), false);
  vector<int> exon_id(trsts.size(), 0);
  int last_start = lpos;
  int last_num = 0;
  long double variance_sum = 0;
  uint64_t variance_count = 0; 
  for(uint64_t i = lpos; i<=rpos; i++){
    bool change = (i == rpos);
    for(int j=0; j<trsts.size(); j++){
      if(in_transcript[j] && i > trsts[j].exons[exon_id[j]].second){
        in_transcript[j] = false;
        change = true; 
        exon_id[j]++;
      }
      else if(!in_transcript[j] && exon_id[j] < trsts[j].exons.size() && i > trsts[j].exons[exon_id[j]].first){
        in_transcript[j] = true;
        change = true; 
      }
    }

    if(change){
      //cout << "Change " << i << " " << last_start << " " << last_num << endl;
      //do calculations from last_start to i-1 if last_num > 0
      if(last_num > 0 && last_start < i){
        int total_read_counts = 0;
        for(int j=last_start; j<i; j++){
          total_read_counts += read_count[j - lpos];
        }
        double avg = (1.0 * total_read_counts)/(i - last_start);
        if(avg == 0) continue;

        double sum = 0;
        for(int j=last_start; j<i; j++){
          sum += abs((read_count[j - lpos] - avg)/avg);
          //sum += pow((read_count[j - lpos] - avg), 2);
          //cout << "j-lpos: " << (j-lpos) << "\tavg: " << avg << "\t" << read_count[j - lpos] << "\t" << (read_count[j - lpos] - avg) << endl;
        }
        //double deviation = pow(sum/(i - last_start), 0.5);
        double variance = sum/(i - last_start);
        variance_count += (i - last_start);
        variance_sum += (i - last_start) * variance;
        //cout << "i:" << i << "\tlast_start:" << last_start << "\tsum:" << sum << "\tvariance_sum: " << variance_sum << "\tvariance_count: " << variance_count << "\t" << (variance_sum/variance_count) << "\t" << variance << endl;
      }
      //cout << "Finished change" << endl;
      last_num = 0;
      for(int j=0; j<trsts.size(); j++) last_num += (in_transcript[j])?1:0;
      last_start = i;
      //cout << "Finished set" << endl;
    }
  } 

  vector<uint64_t> rtn;
	rtn.push_back(mapped);
  rtn.push_back(total_distance);
  rtn.push_back(norm_distance);
  rtn.push_back(mapped_pairs);
  rtn.push_back(reads);
  rtn.push_back(pairs);
  rtn.push_back((variance_count == 0)?9999999999:(int)1000 * variance_sum/variance_count);
      //cout << "Finished loops (" << ((variance_count == 0)?9999999999:(int)variance_sum/variance_count) << ")" << endl;
	return rtn;
}
