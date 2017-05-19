// written and debugged by Joel Graber, Senior Research Associate at the
// Center for Advanced Biotechnology, Boston University
// All rights reserved by The Trustees of Boston University, 1999-2000.
//
//////////////////////////////////////////////////////////////////////////////////////////////
// fastaSeq.h - header file for class FASTAsequence, a class for obtaining and storing
//   FASTA format sequences - initial creation 14 October 1999, Joel Graber
//////////////////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <strstream>
#include <ctype.h>
#include <math.h>
#include <algorithm>
#include "../lib/fastaSeq.h"
#include "../lib/util.h"
#include "../lib/logFile.h"
#include "../lib/filter.h"
typedef vector<unsigned> uVec;
//////////////////////////////////////////////////////////////////////////////////////////////
//
bool ValidSequence::operator()(const string& s) {
  for(unsigned i=0;i<s.length();++i)
    if(!isalpha(s[i]) && !isspace(s[i]) &&(s[i]!='-')) return false;
  // only alphabetic, '-' or white space okay
  return true; // all were okay
}
//////////////////////////////////////////////////////////////////////////////////////////////
bool idxCheck(std::istream& xf,const string& id,int& start,int& N,int& fOffset,string& h) {
  xf.clear(); xf.seekg(0); 
  //  if(xf.fail()) { h=string("ERROR: immediate fail"); return false; }
  //
  bool keyFound=false,multiHits=false;
  static char lineIn[5001]; xf.getline(lineIn,5000);
  unsigned sLength,sLines,lineL; 
  while(!xf.eof()) {
    strVec fields= split()(string(lineIn),'+');
    bool match= (static_cast<int>(fields[0].find(id))>=0);
    //    if(match) WriteToLog()(string("key found: ")+id+string(" in header: ")+fields[0]);
    //    if(match) std::cout << fields[0] << std::endl;
    multiHits= (match&&keyFound) || multiHits;
    keyFound=  match||keyFound;
    if(match && !multiHits) {
      h= fields[0];
      fOffset= atoi(fields[1].c_str()); sLength= atoi(fields[2].c_str());
      sLines=  atoi(fields[3].c_str()); lineL=   atoi(fields[4].c_str());
    }
    xf.getline(lineIn,5000);
  }
  bool noSub= (start<0 && N<0);
  if(start<0) start=0;
  if(N<0) N=sLength-start;
  if(multiHits) { 
    h=string("ERROR:multi-hits"); WriteToLogs()(string("multiple idx file hits for id: ")+id); return false; }
  if(!keyFound) { 
    h=string("ERROR:no hits"); WriteToLogs()(string("no idx file hits for id: ")+id); return false; }
  if(start>sLength) { 
    h=string("ERROR:invalid start point"); WriteIntToErr()(string("start past end of seq: ")+id,start); return false; }
  if((start+N)>sLength) { 
    h=string("ERROR:invalid end point"); WriteIntToErr()(string("stop past end of seq: "),start+N); return false; }
  // now get the true offset to be used
  int addToOffset= start/lineL; // (ends of line up to the start position)
  fOffset+= start+addToOffset+h.length()+2;
  if(!noSub) {
    char buffer[51]; std::strstream buff(buffer,50,std::ios::out);
    buff << ";subsequence: " << start << "->" << (start+N) << '\0';
    h+= string(buffer);
  }
  return true;
} 
//////////////////////////////////////////////////////////////////////////////////////////////
bool GetSeq(std::istream& seqF,int start,int N,string& seq_) {
  seq_.clear(); if(!seqF) return false;
  seqF.seekg(start,std::ios::beg);
  if(seqF.fail() || !seqF) return false;
  int count=0;
  while(count<N && !seqF.eof()) {
    char ch; seqF.get(ch);
    if(ch>='A'&&ch<='Z') { ++count; seq_ += ch; }
    else if (ch>='a'&&ch<='z') { ++count; seq_ += (ch + static_cast<char>('A'-'a')); }
  }
  if(count < N) return false;
  return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
TJfastaSequence::TJfastaSequence(const string& s, const string& h) {
  seq_.reserve(s.length()); head_.reserve(h.length());
  if(ValidSequence()(s)) {
    head_=h;
    for(unsigned i=0;i<s.length();++i)
      if(!isspace(s[i])) {
	if(s[i]=='-')                     seq_ += 'N';
	else if((s[i]>='A')&&(s[i]<='Z')) seq_ += s[i];
	else                              seq_ += (s[i]+'A'-'a');
      }
  } else {
    WriteToLogs()(string("Invalid FASTA sequence: ") + s);
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
TJfastaSequence::TJfastaSequence(const string& fn,const string& in,const string& id,int start,
				 int N):head_(string("not found")),seq_() {
  std::ifstream idxFile(in.c_str(),std::ios::in);
  if(idxFile.fail()||!idxFile){ head_=string("error opening idxFile"); return; }
  int trueStart;
  if(!idxCheck(idxFile,id,start,N,trueStart,head_)) {
    head_ += string(";id= ")+id+string(";idxFile= ")+in; return;
  }
  idxFile.close();
  std::ifstream seqFile(fn.c_str(),std::ios::in);
  if(seqFile.fail()||!seqFile){ 
    head_=string("error opening seqFile"); return; 
  }
  if(!GetSeq(seqFile,trueStart,N,seq_)) {
    head_=string("error reading sequence from ")+fn; return;
  }
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
TJfastaSequence::TJfastaSequence(std::istream& seqF,std::istream& idxF,const string& id,int start,
				 int N):head_(string("not found")),seq_() {
  int trueStart;
  if(!idxCheck(idxF,id,start,N,trueStart,head_)) {
    head_+= string(";id= ")+id; return;
  }
  if(!GetSeq(seqF,trueStart,N,seq_)) {
    head_=string("error reading sequence"); return;
  }
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
string TJfastaSequence::Sequence(unsigned i,unsigned L) const {
  if(i >= seq_.length()) {
    WriteIntToErr()(string("Invalid FASTA sequence index in ")+head_+string(": "),i);
    WriteIntToLog()(string("Invalid FASTA sequence index in ")+head_+string(": "),i);
    return string();
  }
  if(!i && !L) return string(seq_);
  return string(seq_.substr(i,L));
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
string TJfastaSequence::NiceSeq(unsigned L) const {
  string ns;
  if((L<10)||(L>100)) L=60;
  for(unsigned i=0;i<seq_.length();i+=L) {
    for(unsigned j=0;((j<L)&&((i+j)<seq_.length()));j+=10)
      ns += (string(" ")+seq_.substr(i+j,10));
    ns += string("\n");
  }
  return string(ns);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
string TJfastaSequence::NiceDustedSeq(unsigned L,char mask, int wind, int lev, int w) const {
  string ns,ds=Dusted(mask,wind,lev,w);
  if((L<10)||(L>100)) L=60;
  for(unsigned i=0;i<ds.length();i+=L) {
    for(unsigned j=0;((j<L)&&((i+j)<ds.length()));j+=10)
      ns += (string(" ")+ds.substr(i+j,10));
    ns += string("\n");
  }
  return string(ns);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
string TJfastaSequence::Dusted(char mask,int wind,int lev,int w) const {
  SetDustWindow()(wind);
  SetDustLevel()(lev);
  SetDustWord()(w);
  SetDustMask()(mask);
  return string(Dust()(seq_));
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
string TJfastaSequence::DPFiltered(int wind,int word,int th,int cut) const {
  SetDPWindow()(wind);
  SetDPWord()(word);
  SetDPThreshold()(th);
  SetDPCutoff()(cut);
  return string(DPFilter()(seq_));
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
TJfastaSequence& TJfastaSequence::operator=(const TJfastaSequence& o) {
  if(*this==o) return *this;
  head_=o.head_; seq_=o.seq_;
  return *this;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
TJfastaSequence TJfastaSequence::Mutate(const string& mSeq,const string& h) const {
  string s;
  unsigned m=0,p=0;
  while(m<mSeq.length()) 
    if((mSeq[m]=='-')||(mSeq[m]=='.')) {
      if(p>=seq_.length()) {
	WriteToErr()(string("TJfastaSequence::Mutate->mut seq longer than seq."));
	WriteToErr()(string("mSeq: ")+mSeq);
	WriteToErr()(string("seq_: ")+seq_);
	return TJfastaSequence();
      }
      s += seq_[p++]; m++;
    } else if(mSeq[m]=='^') { 
      if(p>=seq_.length()) {
	WriteToErr()(string("TJfastaSequence::Mutate->mut seq longer than seq."));
	WriteToErr()(string("mSeq: ")+mSeq);
	WriteToErr()(string("seq_: ")+seq_);
	return TJfastaSequence();
      }
      p++; m++; 
    } else if(mSeq[m]=='[') {
      if(p>=seq_.length()) {
	WriteToErr()(string("TJfastaSequence::Mutate->mut seq longer than seq."));
	WriteToErr()(string("mSeq: ")+mSeq);
	WriteToErr()(string("seq_: ")+seq_);
	return TJfastaSequence();
      }
      while(mSeq[++m]!=']') s += mSeq[m];
      m++;
    } else {
      s += mSeq[m++]; p++;
    }
  if(p<seq_.length()) {
    WriteToErr()(string("TJfastaSequence::Mutate->mut seq shorter than seq_"));
    WriteToErr()(string("mSeq: ")+mSeq);
    WriteToErr()(string("seq_: ")+seq_);
    return TJfastaSequence(s,h);
  }
  return TJfastaSequence(s,h);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
unsigned TJfastaSequences::LoadSequences(const string& fName,bool mut) {
  std::ifstream seqFile; seqFile.open(fName.c_str(),std::ios::in);//|ios::nocreate);
  if(seqFile.fail()|| !seqFile) {
    WriteToErr()(string("Unable to open file: ")+fName);
    WriteToLog()(string("Unable to open file: ")+fName);
    return 0;
  }
  char *line=new char[maxLineLength+1];
  seqFile.getline(line,maxLineLength);  unsigned count=0,lastStdSeq=0;
  while(!seqFile.eof()) {
    ++count;
    while((line[0]!='>')&& !seqFile.eof()) seqFile.getline(line,maxLineLength);
    if(line[0]=='#' || line[0]==';') continue;
    string h(line),s; s.reserve(500000);
    if(h.length()) {
      h= h.substr(1,h.length()-1);
      while(isspace(h[h.length()-1])) h=h.substr(0,h.length()-1);
      seqFile.getline(line,maxLineLength);
      while((line[0]!='>')&& !seqFile.eof()) {
	s += string(line);
	seqFile.getline(line,maxLineLength);
      }
      if((string(line).length())&&(seqFile.eof())) s += string(line);
      TJfastaSequence faSeq;
      if((h[0]!='+')|| !mut) {
	faSeq=TJfastaSequence(s,h);
	if(mut && faSeq.length()) lastStdSeq= size();
      } else 
	faSeq= (*this)[lastStdSeq].Mutate(s,h); //faSeq= at(lastStdSeq).Mutate(s,h);
      
      if(faSeq.length()) push_back(faSeq);
      else {
	WriteIntToErr()(string("Invalid Sequence ")+h+string(", # "),count);
	WriteIntToLog()(string("Invalid Sequence ")+h+string(", # "),count);
      }
    }
  }
  delete[] line;
  seqFile.close();
  double avg= static_cast<double>(TotalCharacters())/static_cast<double>(size());
  WriteToLog()(string("FASTA File read in: ")+fName);
  WriteIntToLog()(string("sequence count: "),size());
  WriteIntToLog()(string("character count: "),TotalCharacters());
  WriteIntToLog()(string("shortest seq: "),ShortestSeq());
  WriteIntToLog()(string("longest seq: "),LongestSeq());
  WriteFloatToLog()(string("average length: "),avg);
  WriteToLog()(string());
  return size();
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
unsigned TJfastaSequences::LoadSequences(const string& fName,vector<unsigned>& offsets) {
  std::ifstream seqFile; seqFile.open(fName.c_str(),std::ios::in);//|ios::nocreate);
  if(seqFile.fail()|| !seqFile) {
    WriteToErr()(string("Unable to open file: ")+fName);
    WriteToLog()(string("Unable to open file: ")+fName);
    return 0;
  }
  char *line=new char[maxLineLength+1];
  seqFile.getline(line,maxLineLength);  unsigned count=0, file_offset=0; 
  offsets.erase(offsets.begin(),offsets.end());
  while(!seqFile.eof()) {
    count++;
    while((line[0]!='>')&& !seqFile.eof()) { 
      file_offset+= string(line).length()+1;
      seqFile.getline(line,maxLineLength);
    }
    string h(line),s; s.reserve(maxLineLength);
    if(h.length()) {
      offsets.push_back(file_offset);
      h= h.substr(1,h.length()-1);
      if(!isalnum(h[h.length()-1])) h=h.substr(0,h.length()-1);
      seqFile.getline(line,maxLineLength); file_offset += string(line).length()+1;
      while((line[0]!='>')&& !seqFile.eof()) {
	s += string(line);
	seqFile.getline(line,maxLineLength);
	file_offset += string(line).length()+1;
      }
      if((string(line).length())&&(seqFile.eof())) s += string(line);
      TJfastaSequence faSeq(s,h);
      if(faSeq.length()) push_back(faSeq);
      else {
	WriteIntToErr()(string("Invalid Sequence Entry # "),count);
	WriteIntToLog()(string("Invalid Sequence Entry # "),count);
      }
    }
  }
  seqFile.close();
  double avg= static_cast<double>(TotalCharacters())/static_cast<double>(size());
  WriteToLog()(string("FASTA File read in: ")+fName);
  WriteIntToLog()(string("sequence count: "),size());
  WriteIntToLog()(string("character count: "),TotalCharacters());
  WriteIntToLog()(string("shortest seq: "),ShortestSeq());
  WriteIntToLog()(string("longest seq: "),LongestSeq());
  WriteFloatToLog()(string("average length: "),avg);
  WriteToLog()(string());
  return size();
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
unsigned TJfastaSequences::LoadFromTBLFile(const string& fName,bool mut) {
  std::ifstream seqFile; seqFile.open(fName.c_str(),std::ios::in);//|ios::nocreate);
  if(seqFile.fail()|| !seqFile) {
    WriteToErr()(string("Unable to open file: ")+fName);
    WriteToLog()(string("Unable to open file: ")+fName);
    return 0;
  }
  char line[maxLineLength+1];
  seqFile.getline(line,maxLineLength);  unsigned count=0, lastStdSeq=0;
  while(!seqFile.eof()) {
    ++count; string lineIn= line;
    strVec fields= split()(lineIn);
    string h=fields[0],s=fields[1];
    if(h.length()) {
      TJfastaSequence faSeq;
      if((h[0]!='+')|| !mut) {
	faSeq=TJfastaSequence(s,h);
	if(mut && faSeq.length()) lastStdSeq= size();
      } else 
	faSeq= (*this)[lastStdSeq].Mutate(s,h); //faSeq= at(lastStdSeq).Mutate(s,h);
      
      if(faSeq.length()) 
	push_back(faSeq);
      else {
	WriteIntToErr()(string("Invalid Sequence Entry # "),count);
	WriteIntToLog()(string("Invalid Sequence Entry # "),count);
      }
    }
    seqFile.getline(line,1000000);
  }
  seqFile.close();
  double avg= static_cast<double>(TotalCharacters())/static_cast<double>(size());
  WriteToLog()(string("TBL File read in: ")+fName);
  WriteIntToLog()(string("sequence count: "),size());
  WriteIntToLog()(string("character count: "),TotalCharacters());
  WriteIntToLog()(string("shortest seq: "),ShortestSeq());
  WriteIntToLog()(string("longest seq: "),LongestSeq());
  WriteFloatToLog()(string("average length: "),avg);
  WriteToLog()(string());
  return size();
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
unsigned TJfastaSequences::TotalCharacters() const {
  unsigned tc=0;
  for(unsigned i=0;i<size();++i) tc += (*this)[i].length(); //tc += at(i).length();
  return unsigned(tc);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
unsigned TJfastaSequences::LongestSeq() const {
  unsigned maxL=0;
  for(unsigned i=0;i<size();++i) if((*this)[i].length()>maxL) maxL= (*this)[i].length(); //maxL=at(i).length();
  return unsigned(maxL);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
unsigned TJfastaSequences::ShortestSeq() const {
  unsigned minL=static_cast<unsigned>(1e10);
  for(unsigned i=0;i<size();++i) if((*this)[i].length()<minL) minL=(*this)[i].length();
  return unsigned(minL);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
double TJfastaSequences::AverageLength() const {
  if(empty()) return 0.0;
  double sum=0.0;
  for(const_iterator it=begin();it!=end();++it) sum+= it->length();
  return sum/static_cast<double>(size());
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
double TJfastaSequences::MedianLength() const {
  if(empty()) return 0;
  unsigned N=size();
  if(N==1) return (*this)[0].length();
  std::vector<unsigned> len; len.reserve(N);
  for(const_iterator it=begin();it!=end();++it) len.push_back(it->length());
  std::sort(len.begin(),len.end());
  return (N%2)?len[(N-1)/2]:0.5*static_cast<double>(len[N/2]+len[N/2-1]);
}  
//////////////////////////////////////////////////////////////////////////////////////////////
//
uVec TJfastaSequences::LengthHist(unsigned nBins,unsigned binW) const {
  if(empty()) return uVec(nBins,0);
  if(nBins<=0) {
    WriteIntToErr()(std::string("non-positive nBins in TJfastaSequences::LengthHist()"),nBins);
    return uVec();
  }    
  if(binW<=0) { 
    WriteIntToErr()(std::string("non-positive binW in TJfastaSequences::LengthHist()"),binW);
    return uVec(nBins,0);
  }
  uVec bins(nBins,0);
  for(const_iterator it=begin();it!=end();++it) {
    unsigned idx= (it->length())/binW;
    if(idx>=nBins) idx=nBins-1;
    ++bins[idx];
  }
  return bins;
}  
//////////////////////////////////////////////////////////////////////////////////////////////
//
iMat TJfastaSequences::CountVector(int L,int order,int pc,bool cumulative,seq_t t) const {
  if(!L) L= LongestSeq();
  if(t==Protein) return iMat(); // need to fix this, but not today.
  int nBases= (t==Protein)?20:4;
  int nWords= pow(nBases,order+1);
  iMat c(L,iVec(nWords,pc));
  for(int i=0;i<size();++i) {
    std::string s((*this)[i].Sequence());
    for(int p=0;p<s.length()-order;++p){
      unsigned idx= SeqToIdx()(s.substr(p,order+1));
      if(idx<nWords) ++c[p][idx];
      if(cumulative) for(int j=idx+1;j<nWords;++j) ++c[p][j];
    }
  }
  return c;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
TJfastaSequences TJfastaSequences::ProfileRandomize(int& seed,int order,int profL,int N,int pcWeight) const {
  TJfastaSequences pSeqs; pSeqs.reserve(size());
  bool fixedSeqL= (profL>0);
  if(profL<0) 
    profL= LongestSeq();
  if(N<0)
    N= size();
  if(seed<0) {
    time_t prt; time(&prt); seed= prt;
  } 
  int nSeqs=size();
  iMat counts= CountVector(profL,order,pcWeight,true);
  GSLrng prand(seed); 
  int nWords=pow(4,order+1),nsWords=pow(4,order);
  pSeqs.reserve(N);
  for(int i=0;i<N;++i) {
    std::string h(std::string("profileRandomSequence_")+IntToString()(i+1));
    int stop= (fixedSeqL)?profL:(*this)[i%nSeqs].length();
    std::string seq; seq.reserve(stop);
    // start by getting order-1 bases from the whole length
    int lastIdx;
    if(order) {
      int cIdx= prand.uniform_int(counts[0][nWords-1]);
      int sIdx=0;
      while(cIdx>counts[0][sIdx] && sIdx<nWords) ++sIdx;
      seq += IdxToSeq()(sIdx,order+1);
      lastIdx= 4*(sIdx%nsWords);
    }
    for(int p=order;p<stop-order;++p) 
      if(!order) {
	int cIdx= prand.uniform_int(counts[p][3]);
	int sIdx=0;
	while(cIdx>counts[p][sIdx] && sIdx<4) ++sIdx;
	seq += DNAbases[sIdx];
      }	else {
	iVec subCounts(4,0);
	int offset= (lastIdx)?counts[p][lastIdx-1]:0;
	for(int k=0;k<4;++k) subCounts[k]= counts[p][lastIdx+k]-offset;
	int cIdx= prand.uniform_int(subCounts[3]);
	int sIdx=0;
	while(cIdx>subCounts[sIdx] && sIdx<4) ++sIdx;
	seq += IdxToSeq()(sIdx,1);
	lastIdx= 4*((lastIdx+sIdx)%nsWords);
      }
    pSeqs.push_back(TJfastaSequence(seq,h));
  }
  return pSeqs;
}
//////////////////////////////////////////////////////////////////////////////////////////////
// random circular permutation, return the offests
TJfastaSequences TJfastaSequences::CircularRandomize(int& seed, iVec& offsets) const {
  if(seed<0) {
    time_t crt; time(&crt); seed= crt+1;
  } 
  //std::cout << "circ 1" << std::endl;
  GSLrng rand(seed);
  offsets.clear(); offsets.reserve(size());
  //std::cout << "circ 2" << std::endl;
  TJfastaSequences cSeqs; cSeqs.reserve(size());
  for(int i=0;i<size();++i) {
    int o= rand.uniform_int((*this)[i].length());
    //std::cout << "circ loop, i=" << i << std::endl;
    std::string s= (*this)[i].Sequence();
    std::rotate(s.begin(),s.begin()+o,s.end());
    std::string h= (*this)[i].Header()+std::string(" [rotated by ")+IntToString()(o)+std::string("]");
    cSeqs.push_back(TJfastaSequence(s,h)); offsets.push_back(o);
  }
  //std::cout << "circ done " << std::endl;
  return cSeqs;
}
//////////////////////////////////////////////////////////////////////////////////////////////
// random circular permutation
TJfastaSequences TJfastaSequences::CircularRandomize(int& seed) const {
  // seed to time if set to -1
  if(seed<0) {
    time_t crt; time(&crt); seed= crt+1;
  } 
  GSLrng rand(seed);
  TJfastaSequences cSeqs; cSeqs.reserve(size());
  for(int i=0;i<size();++i) {
    int o= rand.uniform_int((*this)[i].length());
    std::string s= (*this)[i].Sequence();
    std::rotate(s.begin(),s.begin()+o,s.end());
    std::string h= (*this)[i].Header()+std::string(" [rotated by ")+IntToString()(o)+std::string("]");
    cSeqs.push_back(TJfastaSequence(s,h));
  }
  return cSeqs;
}
      
