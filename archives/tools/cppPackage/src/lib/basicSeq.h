// written and debugged by Joel Graber, Senior Research Associate at the
// Center for Advanced Biotechnology, Boston University
// All rights reserved by The Trustees of Boston University, 1999.
//
#ifndef __BASICSEQ_H
#define __BASICSEQ_H
//////////////////////////////////////////////////////////////////////////////////////////////
#ifndef _STRING_
#include <string>
#endif
#ifndef _VECTOR_
#include <vector>
#endif
#ifndef _MAP_
#include <map>
#endif
#ifndef _SET_
#include <set>
#endif
#ifndef _FUNCTIONAL_
#include <functional>
#endif
using std::string;
using std::vector;
using std::map;
using std::less;
using std::set;
const string DNAbases("CTAGRYWSMKBDHVN"),rcDNAbases("GATCYRWSKMVHDBN");
const string RNAbases("CUAGRYWSMKBDHVN"),rcRNAbases("GAUCYRWSKMVHDBN");
const string Allbases("CTUAGRYWSMKBDHVN");
const string ProteinBases("ARNDCQEGHILKMFPSTWYVBZX");
const int DNAindex[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //0-15
			-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //16-31
			-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //32-47
			-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //48-63
			//  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O
			-1, 2,10, 0,11,-1,-1, 3,12,-1,-1, 9,-1, 8,14,-1, //64-79
			//  Q  R  S  T  U  V  W  X  Y  Z 
			-1,-1, 4, 7, 1,-1,13, 6,-1, 5,-1,-1,-1,-1,-1,-1, //80-95
			//  a  b  c  d  e  f  g  h  i  j  k  l  m  n  o
			-1, 2,10, 0,11,-1,-1, 3,12,-1,-1, 9,-1, 8,14,-1, //96-111
			//  q  r  s  t  u  v  w  x  y  z
			-1,-1, 4, 7, 1,-1,13, 6,-1, 5,-1,-1,-1,-1,-1,-1 };
const int RNAindex[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //0-15
			-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //16-31
			-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //32-47
			-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //48-63
			//  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O
			-1, 2,10, 0,11,-1,-1, 3,12,-1,-1, 9,-1, 8,14,-1, //64-79
			//  Q  R  S  T  U  V  W  X  Y  Z 
			-1,-1, 4, 7,-1, 1,13, 6,-1, 5,-1,-1,-1,-1,-1,-1, //80-95
			//  a  b  c  d  e  f  g  h  i  j  k  l  m  n  o
			-1, 2,10, 0,11,-1,-1, 3,12,-1,-1, 9,-1, 8,14,-1, //96-111
			//  q  r  s  t  u  v  w  x  y  z
			-1,-1, 4, 7,-1, 1,13, 6,-1, 5,-1,-1,-1,-1,-1,-1 };
const int Protindex[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //0-15
			 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //16-31
			 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //32-47
			 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //48-63
			 //  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O
			 -1, 0,20, 4, 3, 6,13, 7, 8, 9,-1,11,10,12, 2,-1, //64-79
			 //  Q  R  S  T  U  V  W  X  Y  Z 
			 14, 5, 1,15,16,-1,19,17,22,18,21,-1,-1,-1,-1,-1, //80-95
			 //  a  b  c  d  e  f  g  h  i  j  k  l  m  n  o
			 -1, 0,20, 4, 3, 6,13, 7, 8, 9,-1,11,10,12, 2,-1, //96-111
			 //  q  r  s  t  u  v  w  x  y  z
			 14, 5, 1,15,16,-1,19,17,22,18,21,-1,-1,-1,-1,-1 };
const unsigned naBaseN=4,prBaseN=20;
typedef map<char,string,std::less<char> > tjTrBaseMap;
typedef vector<tjTrBaseMap> tjTrBaseMapVector;
typedef map<string,string,std::less<string> > tjCodonMap;
typedef vector<string> strVec;
typedef vector<strVec> strMat;
typedef set<string,less<string> > strSet;
typedef vector<char> cVec;
typedef vector<cVec> cMat;
typedef vector<bool> bVec;
typedef vector<int> iVec;
typedef vector<iVec> iMat;
typedef vector<double> dVec;
typedef enum { DNA, RNA, Protein } seq_t;
typedef map<string,int,std::less<string> > strIntMap;
const int maxKvectorK= 16;
//////////////////////////////////////////////////////////////////////////////////////////////
class  toUpperNoWS:public std::binary_function<string,bool,string> {
 public:
  string operator()(const string& s,bool onlyAlpha=false) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class  toLowerNoWS:public std::unary_function<string,string> {
 public:
  string operator()(const string& s) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class  ValidDNASequence:public std::binary_function<string,bool,bool> {
 public:
  bool operator()(const string& s,bool wc=false) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class  ValidRNASequence:public std::binary_function<string,bool,bool> {
 public:
  bool operator()(const string& s,bool wc=false) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class  PurgeDNASequence:public std::binary_function<string,bool,string> {
 public:
  string operator()(const string& s,bool wc=false) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class  PurgeRNASequence:public std::binary_function<string,bool,string> {
 public:
  string operator()(const string& s,bool wc=false) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class  ValidProteinSequence:public std::binary_function<string,bool,bool> {
 public:
  bool operator()(const string& s,bool wc=false) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class  ReverseComplement:public std::binary_function<string,bool,string> {
 public:
  string operator()(const string& s,bool rna=false) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class  SelfComplement:public std::binary_function<bool,string,bool> {
 public:
  bool operator()(const string& s,bool rna=false) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class  Complement:public std::binary_function<string,bool,string> {
 public:
  string operator()(const string& s,bool rna=false) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class  Transcribe:public std::binary_function<string,bool,string> {
 public:
  string operator()(const string& s,bool wc=false) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class  SetTranscrBases:public std::unary_function<tjTrBaseMap,bool> {
 public:
  bool operator()(const tjTrBaseMap& tb) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class  Replicate:public std::binary_function<string,string,string> {
 public:
  string operator()(const string& s,const string& pr=string()) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class  SetReplicateBases:public std::unary_function<tjTrBaseMap,bool> {
 public:
  bool operator()(const tjTrBaseMap& tb) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class SeqToMdx:public std::unary_function<string,int> {
 public:
  int operator()(const string& s) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class SeqToIdx:public std::binary_function<string,bool,int> {
 public:
  int operator()(const string& s,bool seqCheck=true) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class SeqToIdxVec:public std::binary_function<std::string,int,iVec> {
public:
	iVec operator()(const std::string& s,int w) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class DNASeqToIdxVec:public std::binary_function<std::string,int,iVec> {
public:
	iVec operator()(const std::string& s,int w) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class  IdxToSeq:public std::binary_function<int,int,string> {
 public:
  string operator()(int,int) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class ExpandWildCards:public std::binary_function<string,char,strVec> {
 public:
  strVec operator()(const string& w,char na) const; 
};
//////////////////////////////////////////////////////////////////////////////////////////////
class NeighborWords:public std::binary_function<string,int,strSet> {
 public:
  strSet operator()(const string& w,int mm) const; 
};
//////////////////////////////////////////////////////////////////////////////////////////////
class GetSubLWords:public std::binary_function<string,unsigned,strVec> {
 public:
  strVec operator()(const string& w,unsigned l) const; 
};
//////////////////////////////////////////////////////////////////////////////////////////////
class SubWords:public std::binary_function<string,unsigned,strIntMap> {
 public:
  strIntMap operator()(const string& w,unsigned m) const; 
};
//////////////////////////////////////////////////////////////////////////////////////////////
class WCSeqMatch:public std::binary_function<string,string,bool> {
 public:
  bool operator()(const string& q,const string& s) const; 
};
//////////////////////////////////////////////////////////////////////////////////////////////
class WCSeqSearch:public std::binary_function<string,string,iVec> {
 public:
  iVec operator()(const string& q,const string& s) const; 
};
//////////////////////////////////////////////////////////////////////////////////////////////
class kVector:public std::binary_function<string,int,iVec> {
 public:
  iVec operator()(const string& s,int k) const; 
};
//////////////////////////////////////////////////////////////////////////////////////////////
class kPosVector:public std::binary_function<string,int,iMat> {
 public:
  iMat operator()(const string& s,int k) const; 
};
//////////////////////////////////////////////////////////////////////////////////////////////
class wordPeriods:public std::binary_function<string,bool,iVec> {
 public:
  iVec operator()(const string& w,bool princ) const;
}; 
//////////////////////////////////////////////////////////////////////////////////////////////
class TwoWordPeriods:public std::binary_function<strVec,bool,iVec> {
 public:
  iVec operator()(const strVec& w,bool princ) const;
}; 
//////////////////////////////////////////////////////////////////////////////////////////////
class SeqLogoString:public std::binary_function<dVec,int,string> {
 public:
  string operator()(const dVec& p,int N) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class TranslationTable {
 private:
  tjCodonMap tMap_;
  std::string label_,desc_;
  TranslationTable():tMap_(),label_(),desc_() {}
 public:
  TranslationTable(const std::string& fname);
  TranslationTable(const TranslationTable& o):tMap_(o.tMap_),label_(o.label_),desc_(o.desc_) {}
  std::string Translate(const std::string& rna,int frame=3,bool terse=true) const;
  std::string Label() const { return label_; }
  std::string Desc() const { return desc_; }
};
//////////////////////////////////////////////////////////////////////////////////////////////
#endif
