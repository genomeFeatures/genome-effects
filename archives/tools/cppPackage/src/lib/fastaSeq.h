//////////////////////////////////////////////////////////////////////////////////////////////
// fastaSeq.h - header file for class FASTAsequence, a class for obtaining and storing
//   FASTA format sequences - initial creation 14 October 1999, Joel Graber
//////////////////////////////////////////////////////////////////////////////////////////////
// written and debugged by Joel Graber, Senior Research Associate at the
// Center for Advanced Biotechnology, Boston University
// All rights reserved by The Trustees of Boston University, 1999-2000.
//
#ifndef __FASTASEQ_H
#define __FASTASEQ_H
#ifndef _STRING_
#include <string>
#endif
#ifndef _VECTOR_
#include <vector>
#endif
#ifndef _FUNCTIONAL_
#include <functional>
#endif
#include <iostream>
#include "basicSeq.h"
using std::string; using std::vector;
typedef vector<int> iVec;
typedef vector<unsigned> uVec;
typedef vector<iVec> iMat;
const int maxLineLength=1000000;
//////////////////////////////////////////////////////////////////////////////////////////////
class  ValidSequence:public std::unary_function<string,bool> {
 public: bool operator()(const string& s);
};
//////////////////////////////////////////////////////////////////////////////////////////////
class  TJfastaSequence {
 protected:
  string seq_,head_;
 public:
  TJfastaSequence():seq_(),head_() {}
  TJfastaSequence(const TJfastaSequence& s):seq_(s.seq_),head_(s.head_) {}
  TJfastaSequence(const string& fn,const string& in,const string& id,int start= -1,int N= -1);
  TJfastaSequence(std::istream& seqF,std::istream& idxF,const string& id,int start= -1,int N= -1);
  TJfastaSequence(const string& s, const string& h);
  string Sequence(unsigned s=0,unsigned L=0) const;
  string NiceSeq(unsigned L=60) const;
  string NiceDustedSeq(unsigned L=60,char mask='N',int wind=64,int lev=20,int w=3) const;
  string Header() const { return string(head_); }
  unsigned length() const { return seq_.length(); }
  string Dusted(char mask='N',int wind=64,int lev=20,int w=3) const;
  string DPFiltered(int wind=64,int word=3,int th=1,int cut=100) const;
  bool operator==(const TJfastaSequence& o) const 
  { return((seq_==o.seq_)&&(head_==o.head_)); }
  TJfastaSequence& operator=(const TJfastaSequence& o);
  TJfastaSequence Mutate(const string& mSeq,const string& h) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class  TJfastaSequences:public vector<TJfastaSequence> {
 public:
  unsigned LoadSequences(const string& fName,bool mut=false);
  unsigned LoadSequences(const string& fName,vector<unsigned>& offsets);
  unsigned LoadFromTBLFile(const string& fName,bool mut=false);
  unsigned TotalCharacters() const;
  unsigned LongestSeq() const;
  unsigned ShortestSeq() const;
  double   AverageLength() const;
  double   MedianLength() const;
  uVec     LengthHist(unsigned nBins,unsigned binW) const;
  iMat CountVector(int L=0,int order=0,int pc=1,bool cumulative=false,seq_t t=DNA) const;
  TJfastaSequences ProfileRandomize(int& seed,int order=0,int profL= -1,int N= -1,int pcWeight=1) const;
  TJfastaSequences CircularRandomize(int& seed, iVec& offsets) const;
  TJfastaSequences CircularRandomize(int& seed) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
#endif

