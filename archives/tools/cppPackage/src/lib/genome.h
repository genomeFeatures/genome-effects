/////////////////////////////////////////////////////////////////////////////////////////////
#ifndef __genome_h
#define __genome_h
#include "util.h"
#include "fastaSeq.h"
#include <string>
#include <utility>
#include <iostream>
using std::string;
/////////////////////////////////////////////////////////////////////////////////////////////
// constrained to be a *contiguous* segment of genomic sequence.
class genomeSegment {
 protected:
  intRange r_;
  string chr_,org_;
  bool str_;
 public:
  genomeSegment(const string& org,const string& chr,int i0,int i1);
  genomeSegment(const string& org,const string& chr,intRange r,bool str);
  genomeSegment(const genomeSegment& gs):chr_(gs.chr_),org_(gs.org_),str_(gs.str_),r_(gs.r_) {}
  genomeSegment():chr_(),org_(),str_(),r_(std::make_pair(-1,-1)) {}
  //
  //  TJfastaSequence Sequence() const;
  string Sequence() const;
  string Junction5Sequence(int internalN,int externalN) const;
  string Junction3Sequence(int internalN,int externalN) const;
  //
  bool merge(const genomeSegment& gs);
  bool overlap(const genomeSegment& gs,bool forceStrand=false) const;
  bool contains(const genomeSegment& gs,bool forceStrand=false) const;
  //
  bool operator==(const genomeSegment& gs) const { return(r_==gs.r_ && chr_==gs.chr_ && org_==gs.org_ && str_==gs.str_); }
  genomeSegment& operator=(const genomeSegment& gs) {
    if(!(*this==gs)) {
      r_=gs.r_; chr_=gs.chr_; org_=gs.org_; str_=gs.str_;
    }
    return *this;
  }    
  int begin() const { return r_.first; }
  int end() const { return r_.second; }
  string Chromosome() const { return chr_; }
  string Organism() const { return org_; }
  bool Sense() const { return !str_; }
};
typedef vector<genomeSegment> gsVec;
/////////////////////////////////////////////////////////////////////////////////////////////
class TU {
 protected:
  string org_,chr_,tID_,gID_,pID_;
  int start_,stop_,utr5L_,utr3L_,cdsL_,tL_;
  bool str_; irVec exons_,introns_;
 public:
  typedef enum{ ensembl,sgd,tair } datatype;
  TU() {}
  TU(const TU& o):org_(o.org_),chr_(o.chr_),tID_(o.tID_),gID_(o.gID_),pID_(o.pID_),start_(o.start_),stop_(o.stop_),
    utr5L_(o.utr5L_),utr3L_(o.utr3L_),str_(o.str_),exons_(o.exons_),introns_(o.introns_),cdsL_(o.cdsL_),tL_(o.tL_) {}
  TU(const string& org,const string& line,datatype dt=ensembl);
  //
  string Sequence() const;
  string CDS() const { return Sequence().substr(utr5L_,cdsL_); }
  gsVec Exons() const;
  gsVec Introns() const;
  strVec IntronSequences() const;
  strVec ExonSequences() const;
  strVec ExonIntronSequences(int eN,int iN) const;
  strVec IntronExonSequences(int iN,int eN) const;
  string UTR5() const;
  string UTR3() const;
  int length() const;
  string Upstream(int N) const;
  string Downstream(int N) const;
  //
  string Organism() const { return org_; }
  string Chromosome() const { return chr_; }
  string TranscriptID() const { return tID_; }
  string GeneID() const { return gID_; }
  string ProteinID() const { return pID_; }
  bool Strand() const { return str_; }
  char Str() const { return str_?'-':'+'; }
  irVec ExonBoundaries() const { return exons_; }
  irVec IntronBoundaries() const { return introns_; }
  int UTR5length() const { return utr5L_; }
  int UTR3length() const { return utr3L_; }
  int nExons() const { return exons_.size(); }
  int begin() const { return start_; }
  int end() const { return stop_; }
  int CDSLength() const { return cdsL_; }
  int transcriptLength() const { return tL_; }
};
/////////////////////////////////////////////////////////////////////////////////////////////
class Repeat:public genomeSegment {
 protected:
  string rID_,stID_,tID_,detailed_,method_;
 public:
  Repeat() {}
  Repeat(const Repeat& o):genomeSegment(o),rID_(o.rID_),stID_(o.stID_),tID_(o.tID_),detailed_(o.detailed_),method_(o.method_) {}
  Repeat(const string& line);
  Repeat(std::istream& is);
  string RptSubType() const { return string(stID_); }
  string RptType() const { return string(tID_); }
  string Method() const { return string(method_); }
  string DetailedID() const { return string(detailed_); }
  string ID() const { return string(rID_); }
  bool operator==(const Repeat& o) {
    return (r_==o.r_ && chr_==o.chr_ && org_==o.org_ && str_==o.str_ && rID_==o.rID_ && tID_==o.tID_ && stID_==o.stID_ && detailed_==o.detailed_ && method_==o.method_);
  }
  Repeat& operator= (const Repeat& o) { 
    if(!(*this==o)) {
      r_=o.r_; chr_=o.chr_; org_=o.org_; str_=o.str_; rID_=o.rID_; stID_=o.stID_; detailed_=o.detailed_; tID_=o.tID_; method_=o.method_; return *this; 
    }
    return *this;
  }
  void SetOrg(const string& o) { org_=o; }
};
typedef vector<Repeat> rptVector;
/////////////////////////////////////////////////////////////////////////////////////////////
#endif
/////////////////////////////////////////////////////////////////////////////////////////////
