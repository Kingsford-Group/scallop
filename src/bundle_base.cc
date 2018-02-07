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
{
  hits.clear();
}

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

vector<uint64_t> bundle_base::hits_on_transcripts(int bundle_index, unordered_map <string,bool> &pair_mapping){
	int mapped = 0;
  float total_distance = 0;
  float norm_distance = 0;
  int mapped_pairs = 0;
  int pairs = 0;
  int reads = 0;
  unordered_map<string,int> hit_index;
  vector<uint64_t> read_count(rpos - lpos + 1, 0);
  for(int i=0; i<hits.size(); i++){
    string name = hits[i].qname;
    name += ((hits[i].flag & 0x40)?".1":".2");
    assert(!((hits[i].flag & 0x40) && (hits[i].flag & 0x80)));
    assert((hits[i].flag & 0x40) || (hits[i].flag & 0x80));
    bool found = false;
		//if(pair_mapping.find(name) == pair_mapping.end()){
      reads++;
			//pair_mapping[name] = false;
    //}

		for(int j=0; j<trsts.size(); j++){
			if(hits[i].maps_to_transcript(trsts[j])){
        //cerr << name << "\t" << bundle_index << "\t" << trsts[j].seqname << "\t" << trsts[j].gene_id << "\t" << trsts[j].transcript_id << endl;
        if(!found){//pair_mapping[name] == false){
          mapped++;
          found = true;
          total_distance += hits[i].nm;
          norm_distance += 1.0 * hits[i].nm / hits[i].qlen;
					//cout << "Mapping " << name ;
          for(int l=0; l<=hits[i].spos.size(); l++){
            uint64_t start = (l==0)?hits[i].pos:low32(hits[i].spos[l-1]);
            uint64_t end = (l==hits[i].spos.size())?hits[i].rpos:high32(hits[i].spos[l]);
						//cout << "\t(" << start << "," << end << ")";
            for(uint64_t k=start; k<=end; k++){
              read_count[k - lpos]++;
            }
          }
					//cout << endl;

        }
        //pair_mapping[name] = true;
      }
		}

    if(hit_index.find(hits[i].qname) != hit_index.end()){
      if(true){//pair_mapping.find(hits[i].qname) == pair_mapping.end()){
        pairs++;
      	//pair_mapping[hits[i].qname] = false;
			}
      int ip = hit_index[hits[i].qname];
      bool go = true;
      for(int j=0; j<trsts.size() && go; j++){
        if(hits[i].maps_to_transcript(trsts[j]) && hits[ip].maps_to_transcript(trsts[j])){
          if(true){//pair_mapping[hits[i].qname] == false){
            mapped_pairs++;
          }
          //pair_mapping[hits[i].qname] = true;
          go = false;
          break;
        }
      }
    }
    hit_index[hits[i].qname] = i;

  }

	/*if(trsts.size() > 0 && trsts[0].transcript_id == "gene.0.0.1"){
		for(int j=0; j<trsts.size(); j++){
			cout << "Transcript " << j << ":";
			for(int exon = 0; exon < trsts[j].exons.size(); exon++){
				cout << " (" << trsts[j].exons[exon].first << "," << trsts[j].exons[exon].second << ")";
			}
			cout << endl;
		}
		for(int j=0; j<hits.size(); j++){
			cout << "Hit " << j << "(" << hits[j].qname << "):";
			for(int l=0; l<=hits[j].spos.size(); l++){
				uint64_t start = (l==0)?hits[j].pos:low32(hits[j].spos[l-1]);
				uint64_t end = (l==hits[j].spos.size())?hits[j].rpos:high32(hits[j].spos[l]);
				cout << " (" << start << "," << end << ")";
			}
			cout << "\t";
			for(int k=0; k<trsts.size(); k++)
			 cout << hits[j].maps_to_transcript(trsts[k]) << " ";
			cout << endl;
		}
		exit(0);
	}*/

  double coverage_sum = 0;
  for(int j=0; j<trsts.size(); j++) coverage_sum += trsts[j].coverage;
  for(int j=0; j<trsts.size(); j++) trsts[j].covratio = trsts[j].coverage/coverage_sum;
  //for(int j=0; j<trsts.size(); j++) cout << trsts[j].covratio << "\t" << trsts[j].coverage << "\t" << coverage_sum << endl;

  vector<bool> in_transcript(trsts.size(), false);
  vector<int> exon_id(trsts.size(), 0);

  int overall_sum = 0;
  int overall_count = 0;
  for(int j=0; j<trsts.size(); j++) in_transcript[j] = (lpos >= trsts[j].exons[0].first);
  for(uint64_t i = lpos; i<=rpos; i++){
    bool in_one = false;
    for(int j=0; j<trsts.size(); j++){
      if(in_transcript[j] && i > trsts[j].exons[exon_id[j]].second){
        in_transcript[j] = false;
        exon_id[j]++;
      }else if(in_transcript[j]){
        in_one = true;
      }else if(!in_transcript[j] && exon_id[j] < trsts[j].exons.size() && i >= trsts[j].exons[exon_id[j]].first){
        in_transcript[j] = true;
        in_one = true;
      }
    }
    if(in_one){
      overall_sum += read_count[i - lpos];
      //cout << "In: " << read_count[i - lpos] << endl;
      overall_count ++;
    }else{
      //cout << "Out:" << read_count[i - lpos] << endl;
    }
  }
  double overall_average = (1.0 * overall_sum) / overall_count;
  int last_start = lpos;
  int last_num = 0;
  long double variance_sum = 0;
	long double variance_predicted_sum = 0;
  uint64_t variance_count = 0;
	double graph_sum = 0;
	for(int j=0; j<trsts.size(); j++){
    exon_id[j] = 0;
    in_transcript[j] = (lpos >= trsts[j].exons[0].first);
    last_num += (in_transcript[j])?1:0;
    //graph_sum += ((in_transcript[j])?1:0) * trsts[j].coverage;
    graph_sum += ((in_transcript[j])?1:0) * trsts[j].covratio * overall_average;
  }
  for(uint64_t i = lpos; i<=rpos; i++){
    bool change = (i == rpos);
    for(int j=0; j<trsts.size(); j++){
			//cout << "i:" << i << "\tj:" << j << "\tin_transcript[j]:" << in_transcript[j] << "\ttrsts[j].exons[exon_id[j]].second:" << trsts[j].exons[exon_id[j]].second << endl;
      if(in_transcript[j] && i > trsts[j].exons[exon_id[j]].second){
        in_transcript[j] = false;
        change = true;
        exon_id[j]++;
      }
      else if(!in_transcript[j] && exon_id[j] < trsts[j].exons.size() && i >= trsts[j].exons[exon_id[j]].first){
        in_transcript[j] = true;
        change = true;
      }
    }

    if(change){
      //do calculations from last_start to i-1 if last_num > 0
      if(last_num > 0 && last_start < i){
        int total_read_counts = 0;
        for(int j=last_start; j<i; j++){
          total_read_counts += read_count[j - lpos];
        }

        double avg = (1.0 * total_read_counts)/(i - last_start);
        if(avg > 0){
	        double sum = 0;
	        for(int j=last_start; j<i-1; j++){
	          variance_sum += abs((read_count[j - lpos] - avg)/avg);
						if(graph_sum != 0)
							variance_predicted_sum += abs((read_count[j - lpos] - graph_sum)/graph_sum);
	        }
	        variance_count += (i - last_start);
				}else if(false){
          cout << "last_start:" << last_start << "\ti:" << i << endl;
          for(int j=0; j<trsts.size(); j++){
            cout << "Transcript " << j << ":";
            for(int exon = 0; exon < trsts[j].exons.size(); exon++){
              cout << " (" << trsts[j].exons[exon].first << "," << trsts[j].exons[exon].second << ")";
            }
            cout << endl;
          }

          for(int j=0; j<hits.size(); j++){
            cout << "Hit " << j << ":";
            for(int l=0; l<=hits[j].spos.size(); l++){
              uint64_t start = (l==0)?hits[j].pos:low32(hits[j].spos[l-1]);
              uint64_t end = (l==hits[j].spos.size())?hits[j].rpos:high32(hits[j].spos[l]);
              cout << " (" << start << "," << end << ")";
            }
            cout << endl;
          }

        }
      }
      //cout << "Finished change" << endl;
      last_num = 0;
      for(int j=0; j<trsts.size(); j++) last_num += (in_transcript[j])?1:0;
			graph_sum = 0;
			//for(int j=0; j<trsts.size(); j++) graph_sum += ((in_transcript[j])?1:0) * trsts[j].coverage;
			for(int j=0; j<trsts.size(); j++) graph_sum += ((in_transcript[j])?1:0) * trsts[j].covratio * overall_average;
      last_start = i;
      //cout << "Finished set" << endl;

      //cout << "graph_sum:" << graph_sum ;
      //for(int j=0; j<trsts.size(); j++) cout << "\t(" << trsts[j].covratio << "," << overall_average << ")";
      //cout << endl;
    }


  }

  //if(trsts.size() > 0 && trsts[0].transcript_id == "gene.0.0.1"){ exit(0); }

  vector<uint64_t> rtn;
	rtn.push_back(mapped);
  rtn.push_back(total_distance);
  rtn.push_back(norm_distance);
  rtn.push_back(mapped_pairs);
  rtn.push_back(reads);
  rtn.push_back(pairs);
  if(trsts.size()>0){
    //cout << "rsts[0].transcript_id:" << trsts[0].transcript_id << "\t" << "variance_sum/variance_count:" << (variance_sum/variance_count) << "\tvariance_sum:" << variance_sum << "\tvariance_count:" << variance_count << endl;
  }
  rtn.push_back((variance_count == 0)?9999999999:(int)1000 * variance_sum/variance_count);
  rtn.push_back((variance_count == 0)?9999999999:(int)1000 * variance_predicted_sum/variance_count);
      //cout << "Finished loops (" << ((variance_count == 0)?9999999999:(int)variance_sum/variance_count) << ")" << endl;
  hit_index.clear();
  trsts.clear();
  return rtn;
}
