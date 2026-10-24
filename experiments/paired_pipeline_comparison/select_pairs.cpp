// Select prespecified paired record IDs; never inspect sequence content to select.
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <array>
using namespace std;
int main(int argc,char**argv) { try {
  if(argc!=10) throw runtime_error("usage: ids r1 r2 out-prefix accession first last deep-mode expected-selected");
  ifstream ids(argv[1]),r1(argv[2]),r2(argv[3]);
  if(!ids||!r1||!r2) throw runtime_error("input unavailable");
  const string prefix=argv[4], accession=argv[5];
  const auto first=stoull(argv[6]),last=stoull(argv[7]);
  const bool deep=stoi(argv[8]); const auto expected=stoull(argv[9]);
  ofstream s1(prefix+".shallow_1.fastq"),s2(prefix+".shallow_2.fastq"),d1,d2;
  if(deep){d1.open(prefix+".deep_1.fastq");d2.open(prefix+".deep_2.fastq");}
  unsigned long long wanted=0,selected=0,shallow=0; int is_shallow=0;
  bool has=bool(ids>>wanted>>is_shallow);
  if(has && (wanted<first||wanted>last)) throw runtime_error("ID outside chunk");
  for(auto i=first;i<=last;i++) {
    array<string,4>a,b;
    for(int j=0;j<4;j++) if(!getline(r1,a[j])||!getline(r2,b[j])) throw runtime_error("truncated pair");
    if(a[0]!="@"+accession+"."+to_string(i)+"/1"||b[0]!="@"+accession+"."+to_string(i)+"/2") throw runtime_error("unexpected pair ID");
    for(auto*rec:{&a,&b}) if((*rec)[1].empty()||(*rec)[1].size()!=(*rec)[3].size()||(*rec)[2].empty()||(*rec)[2][0]!='+') throw runtime_error("invalid FASTQ record");
    if(has&&i==wanted) {
      if(deep) for(int j=0;j<4;j++){d1<<a[j]<<'\n';d2<<b[j]<<'\n';}
      if(is_shallow) {for(int j=0;j<4;j++){s1<<a[j]<<'\n';s2<<b[j]<<'\n';} shallow++;}
      selected++; auto previous=wanted; has=bool(ids>>wanted>>is_shallow);
      if(has && (wanted<=previous||wanted>last)) throw runtime_error("unordered or out-of-range selection");
    }
  }
  string extra;
  if(getline(r1,extra)||getline(r2,extra)||has||selected!=expected) throw runtime_error("record count mismatch");
  s1.close();s2.close(); if(deep){d1.close();d2.close();}
  if(!s1||!s2||(deep&&(!d1||!d2))) throw runtime_error("output write failure");
  cout<<"{\"selected\":"<<selected<<",\"shallow\":"<<shallow<<",\"scanned\":"<<(last-first+1)<<"}\n";
  return 0;
}catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
