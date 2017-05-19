// written and debugged by Joel Graber, Senior Research Associate at the
// Center for Advanced Biotechnology, Boston University
// All rights reserved by The Trustees of Boston University, 1999.
//
#include <ctype.h>
#include <math.h>
#ifdef INTEL
#include <mathimf.h>
#endif
#include <strstream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "../lib/basicSeq.h"
#include "../lib/logFile.h"
#include "../lib/util.h"
//////////////////////////////////////////////////////////////////////////////////////////////
// global variable to hold the standard transcription units;
tjTrBaseMap TranscrBases,ReplicBases;
//////////////////////////////////////////////////////////////////////////////////////////////
//
std::string toUpperNoWS::operator()(const std::string& s,bool onlyAlpha) const {
  std::string uc; uc.reserve(s.length());
  char offset= 'A'-'a';
  for(unsigned i=0;i<s.length();++i)
    if((s[i]>='a')&&(s[i]<='z'))               
      uc += (s[i]+offset);
    else if((!onlyAlpha && !isspace(s[i])) || ((s[i]>='A')&&(s[i]<='Z')))      
      uc += s[i];
  return uc;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
std::string toLowerNoWS::operator()(const std::string& s) const {
  std::string uc; uc.reserve(s.length());
  char offset= 'a'-'A';
  for(unsigned i=0;i<s.length();++i)
    if((s[i]>='A')&&(s[i]<='Z')) uc += (s[i]+offset);
    else if(!isspace(s[i]))      uc += s[i];
  return uc;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
bool ValidDNASequence::operator()(const std::string& s,bool wc) const {
  std::string ss=toUpperNoWS()(s);
  for(unsigned i=0;i<ss.length();++i)
    if(ss[i]=='-') ss[i] = 'N'; // allow for gaps
    else {
      unsigned p= DNAbases.find(ss[i]);
      if      ( (p>=naBaseN) && !wc )  return false;
      else if   (p>=DNAbases.length()) return false;
    }
  return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
bool ValidRNASequence::operator()(const std::string& s,bool wc) const {
  std::string ss=toUpperNoWS()(s);
  for(unsigned i=0;i<ss.length();++i)
    if(s[i]=='-') ss = 'N';
    else {
      unsigned p= RNAbases.find(ss[i]);
      if     ( (p>=naBaseN) && !wc )  return false;
      else if  (p>=RNAbases.length()) return false;
    }
  return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
std::string PurgeDNASequence::operator()(const std::string& s,bool wc) const {
  std::string ss=toUpperNoWS()(s,true);
  std::string ps; ps.reserve(ss.length());
  for(unsigned i=0;i<ss.length();++i) {
    unsigned p= DNAbases.find(ss[i]);
    if( (p<naBaseN) || (wc && p<DNAbases.length()) ) ps += ss[i];
  }
  return ps;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
std::string PurgeRNASequence::operator()(const std::string& s,bool wc) const {
  std::string ss=toUpperNoWS()(s,true);
  std::string ps; ps.reserve(ss.length());
  for(unsigned i=0;i<ss.length();++i) {
    unsigned p= RNAbases.find(ss[i]);
    if( (p<naBaseN) || (wc && p<RNAbases.length()) ) ps += ss[i];
  }
  return ps;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
bool ValidProteinSequence::operator()(const std::string& s,bool wc) const {
  std::string ss=toUpperNoWS()(s);
  for(unsigned i=0;i<ss.length();++i) {
    if(ss[i]=='-') ss[i]='X';
    else {
      unsigned p=ProteinBases.find(ss[i]);
      if    ( (p>=prBaseN) &&  !wc )     return false;
      else if (p>=ProteinBases.length()) return false;
    }
  }
  return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
std::string ReverseComplement::operator()(const std::string& s,bool rna) const {
  if(!rna&&!ValidDNASequence()(s,true)) return std::string();
  else if(rna&&!ValidRNASequence()(s,true)) return std::string();
  std::string ss=toUpperNoWS()(s);
  std::string rc; rc.reserve(ss.length());
  for(unsigned i=ss.length();i>0;--i)
    if(!rna) rc += rcDNAbases[DNAbases.find(ss[i-1])];
    else     rc += rcRNAbases[RNAbases.find(ss[i-1])];
  return rc;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
bool SelfComplement::operator()(const std::string& s,bool rna) const {
  if(!rna&&!ValidDNASequence()(s,true)) return false;
  else if(rna&&!ValidRNASequence()(s,true)) return false;
  if(s.length() % 2) return false; // odd length sequences can't be self-complementary
  std::string ss=toUpperNoWS()(s); bool sc=true;
  for(unsigned i=ss.length();(i>ss.length()/2)&&sc;--i)
    if(!rna) sc = rcDNAbases[DNAbases.find(ss[i-1])]==ss[ss.length()-i];
    else     sc = rcRNAbases[RNAbases.find(ss[i-1])]==ss[ss.length()-i];
  return sc;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
std::string Complement::operator()(const std::string& s,bool rna) const {
  if(!rna&&!ValidDNASequence()(s,true)) return std::string();
  else if(rna&&!ValidRNASequence()(s,true)) return std::string();
  std::string ss=toUpperNoWS()(s);
  std::string rc; rc.reserve(ss.length());
  for(unsigned i=0;i<ss.length();++i)
    if(!rna) rc += rcDNAbases[DNAbases.find(ss[i])];
    else     rc += rcRNAbases[RNAbases.find(ss[i])];
  return rc;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
std::string Transcribe::operator()(const std::string& s,bool wc) const {
  if(!ValidDNASequence()(s,wc)) return std::string();
  if(!TranscrBases.size()) {  // if the transcription bases haven't been set, set them to default
    TranscrBases['C']=std::string("rC"); TranscrBases['T']=std::string("rU");
    TranscrBases['A']=std::string("rA"); TranscrBases['G']=std::string("rG");
  }
  std::string ss=toUpperNoWS()(s);
  std::string tr; tr.reserve(ss.length()*3/2);
  tr += std::string("5ppp-");  // tri-phosphate at 5' most end
  for(unsigned i=0;i<ss.length();++i) {
    unsigned baseIdx= DNAbases.find(ss[i]);
    if (baseIdx <= 3) 
      tr += TranscrBases[ss[i]];
    else if (baseIdx < DNAbases.length()) 
      tr += ss[i];
    else {
      WriteToErr()(std::string("Illegal DNA base in Transcribe."));
      WriteToLog()(std::string("Illegal DNA base in Transcribe."));
      return std::string(); // only basic and defined DNA
    }
  }
  return tr;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
bool SetTranscrBases::operator()(const tjTrBaseMap& tb) const {
  if(tb.size()!=4) return false;
  TranscrBases= tb;
  return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
std::string Replicate::operator()(const std::string& s,const string& pr) const {
  if(!ValidDNASequence()(s,true)) return std::string();
  if(!ReplicBases.size()) {  // if the replication bases haven't been set, set them to default
    ReplicBases['C']=std::string("dC"); ReplicBases['T']=std::string("dU");
    ReplicBases['A']=std::string("dA"); ReplicBases['G']=std::string("dG");
  }
  std::string ss=toUpperNoWS()(s);
  std::string tr; tr.reserve(ss.length()*3/2);
  if(pr.length())	tr = pr;
  else			tr += std::string("5dp-");  // dephosphorylation at 5' most end
  for(unsigned i=pr.length();i<ss.length();++i) {
    unsigned baseIdx= DNAbases.find(ss[i]);
    if (baseIdx <= 3) 
      tr += ReplicBases[ss[i]];
    else if (baseIdx < DNAbases.length()) 
      tr += ss[i];
    else {
      WriteToErr()(std::string("Illegal DNA base in Replicate."));
      WriteToLog()(std::string("Illegal DNA base in Replicate."));
      return std::string(); // only basic and defined DNA
    }
  }
  return tr;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
bool SetReplicateBases::operator()(const tjTrBaseMap& tb) const {
  if(tb.size()!=4) return false;
  ReplicBases= tb;
  return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
int SeqToMdx::operator()(const std::string& s) const {
  if(s.length()>15) return -1;
  int idx=0,rdx=0,L=s.length(); 
  for(unsigned i=0;i<L;++i) {
    if((unsigned)DNAbases.find(s[i])>3) return -1;
    else idx= 4*idx + DNAbases.find(s[i]);
    if((unsigned)rcDNAbases.find(s[L-i-1])>3) return -1;
    else rdx= 4*rdx + rcDNAbases.find(s[L-i-1]);
  }
  return (idx<rdx)?idx:rdx;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
int SeqToIdx::operator()(const std::string& s,bool checkSeq) const {
  if(s.length()>15) return -1;
  bool isDNA=true;
  if(checkSeq) {
    isDNA= ValidDNASequence()(s,false);
    bool isRNA= !DNA && ValidRNASequence()(s,false);
    if(!isDNA && !isRNA) return -1;
  }
  const std::string& bases= isDNA?DNAbases:RNAbases;
  int idx=0;
  for(unsigned i=0;i<s.length();++i)
    idx= 4*idx + bases.find(s[i]);
  return idx;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
iVec SeqToIdxVec::operator()(const std::string& s,int n) const {
  int L=s.length(); iVec iv(L-n+1,-1);
  if(n>15 || n<1) return iv;
  if(n>L) return iv; 
  bool isDNA= ValidDNASequence()(s,true); // allow wildcards somewhere in the sequence
  bool isRNA= !DNA && ValidRNASequence()(s,true);
  if(!isDNA && !isRNA) return iv;
  const std::string& bases= isDNA?DNAbases:RNAbases;
  unsigned max=4;
  for(int i=1;i<n;++i) max *=4;
  //unsigned max=pow(static_cast<double>(4),static_cast<double>(n));
  int mask=max-1,Lmn=L-n;
  //get the first word
  int start=0;
  int idx= SeqToIdx()(s.substr(start,n));
  if(iv.size()==1) iv[0]=idx;
  while(idx<0 && start<=Lmn) { 
    idx= SeqToIdx()(s.substr(++start,n));
  }	
  if(start<Lmn) iv[start]= idx;
  while(++start <= Lmn) {
    idx <<= 2;
    idx &= mask;
    int bIdx= bases.find(s[start+n-1]);
    if(bIdx>=0 && bIdx<4) 
      idx += bIdx;
    else {
      start += n;
      idx = -1;
      while(idx < 0 && start <= Lmn) {
	idx = SeqToIdx()(s.substr(start,n));
	if(idx < 0) ++start;
      }
    }
    if(start <= Lmn && idx>=0) iv[start]= idx;
  }
  return iv;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
iVec DNASeqToIdxVec::operator()(const std::string& s,int n) const {
  int L=s.length(); iVec iv(L-n+1,-1);
  if(n>15 || n<1) return iv;
  if(n>L) return iv; 
  unsigned max=4;
  for(int i=1;i<n;++i) max *=4;
  //unsigned max=pow(static_cast<double>(4),static_cast<double>(n));
  int mask=max-1,Lmn=L-n;
  //get the first word
  int start=0;
  int idx= SeqToIdx()(s.substr(start,n));
  if(iv.size()==1) iv[0]=idx;
  while(idx<0 && start<=Lmn) { 
    idx= SeqToIdx()(s.substr(++start,n));
  }	
  if(start<Lmn) iv[start]= idx;
  while(++start <= Lmn) {
    idx <<= 2;
    idx &= mask;
    int bIdx= DNAbases.find(s[start+n-1]);
    if(bIdx>=0 && bIdx<4) 
      idx += bIdx;
    else {
      start += n;
      idx = -1;
      while(idx < 0 && start <= Lmn) {
	idx = SeqToIdx()(s.substr(start,n),false);
	if(idx < 0) ++start;
      }
    }
    iv[start]= idx;
  }
  return iv;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
std::string IdxToSeq::operator()(int idx,int N) const {
  std::string s; s.reserve(N);
  if(idx<0) 
    for(int i=0;i<N;++i) 
      s += 'N';
  else {
    int div=(N==1)?1:static_cast<int>(pow(4,N-1));
    while(N>0) {
      s += DNAbases[idx/div];
      idx = idx% div;
      div /= 4; --N;
    }
  }
  return s;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
strVec ExpandWildCards::operator()(const std::string& w,char na) const {
  if((na>='a')&&(na<='z')) na += 'A'-'a';
  if((na!='D')&&(na!='R')) {
    std::string t; t+= na;
    WriteToErr()(std::string("Invalid nucleic acid designation in ExpandWildCards: ")+t); 
    return strVec();
  }
  if((na=='D')&&!ValidDNASequence()(w,true)) {
    WriteToErr()(std::string("Invalid DNA sequence in ExpandWildCards: ")+w);
    return strVec();
  }	
  if((na=='R')&&!ValidDNASequence()(w,true)) {
    WriteToErr()(std::string("Invalid RNA sequence in ExpandWildCards: ")+w);
    return strVec();
  }
  cMat varBases; varBases.reserve(w.length());
  char ut='T'; if(na=='R') ut='U'; unsigned count=1;
  for(unsigned i=0;i<w.length();++i) {
    cVec temp;
    switch(w[i]) {
    case 'C': case 'T': case 'A': case 'G': case 'U': temp.push_back(w[i]); break;
      // purine-pyrimidine
    case 'R': temp.push_back('A'); temp.push_back('G'); break;
    case 'Y': temp.push_back('C'); temp.push_back(ut); break;
      // meto-keno
    case 'M': temp.push_back('C'); temp.push_back('A'); break;
    case 'K': temp.push_back(ut); temp.push_back('G'); break;
      // weak-strong
    case 'W': temp.push_back(ut); temp.push_back('A'); break;
    case 'S': temp.push_back('C'); temp.push_back('G'); break;
      // three-letters
    case 'B': temp.push_back('C'); temp.push_back(ut); temp.push_back('G'); break;
    case 'D': temp.push_back(ut); temp.push_back('A'); temp.push_back('G'); break;
    case 'H': temp.push_back('C'); temp.push_back(ut); temp.push_back('A'); break;
    case 'V': temp.push_back('C'); temp.push_back('A'); temp.push_back('G'); break;
      // all-letters
    case 'N': temp.push_back('C'); temp.push_back(ut); temp.push_back('A'); temp.push_back('G'); break;
    }
    varBases.push_back(temp); count *= temp.size();
  }
  strVec expWords; expWords.reserve(count); expWords.push_back(std::string());
  for(unsigned j=0;j<varBases.size();++j) {
    strVec tExpWords=expWords; expWords.erase(expWords.begin(),expWords.end());
    for(unsigned k=0;k<tExpWords.size();++k)
      for(unsigned m=0;m<varBases[j].size();++m) {
	std::string word(tExpWords[k]); 
	word += varBases[j][m];
	expWords.push_back(word);
      }
  }
  return strVec(expWords);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
strSet NeighborWords::operator()(const std::string& w,int mm) const {
  if((mm<0)||(mm>w.length())) {
    WriteIntToErr()(std::string("invalid mm in NeighborWords:"),mm);
    return strSet();
  }
  strSet nWords; nWords.insert(w);
  for(unsigned i=1;i<=mm;++i){
    std::vector<bool> wc; 
    for(unsigned j=0;j<w.length();++j) 
      if(j<i) wc.push_back(true);
      else    wc.push_back(false);
    bool cont=true;
    while(cont) {
      std::string tw(w); for(unsigned j=0;j<w.length();++j) if(wc[j]) tw[j]='N';
      strVec nw= ExpandWildCards()(tw,'D');
      for(unsigned j=0;j<nw.size();++j) nWords.insert(nw[j]);
      // now increment the boolean vector
      char k=wc.size()-1; bool foundFalse=false,Found=false; unsigned char postCount=0;
      while(!Found&&(k>=0)) {
	if(!foundFalse && !wc[k]) foundFalse=true;
	else if(!foundFalse) ++postCount;
	else if(wc[k]) Found=true;
	if(!Found) --k;
      }
      if(Found) {
	wc[k]=false; wc[k+1]=true;
	for(int l=k+2;l<wc.size();++l)
	  if((l-k-2)<postCount) wc[l]=true;
	  else                  wc[l]=false;
      } else
	cont= false;
    }
  }
  return nWords;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
strVec GetSubLWords::operator()(const std::string& w,unsigned l) const {
  if((l>w.length())||!l) {
    WriteIntToErr()(std::string("Invalid l in GetSubLWords: "),l);
    return strVec();
  }
  if(!w.length()) {
    WriteToErr()(std::string("Zero-Length word in GetSubLWords: "));
    return strVec();
  }
  strVec subLWords;
  if(l==1) {
    for(unsigned i=0;i<w.length();++i) {
      std::string t; t+= w[i];
      subLWords.push_back(t);
    }
  } else if(l==w.length()) {
    subLWords.push_back(w);
  } else {
    unsigned mn=l,mx=w.length()-l;
    if(mx<mn) { mn=w.length()-l; mx=l; }
    unsigned count=1;
    for(unsigned i0=mx;i0<w.length();++i0) count *= (i0+1);
    for(unsigned j0=mn;j0>1;--j0) count /= j0;
    bVec onOff; 
    for(unsigned i1=0;i1<l;++i1) onOff.push_back(true);
    for(unsigned j1=l;j1<w.length();++j1) onOff.push_back(false);
    for(unsigned k=0;k<count;++k) {
      std::string lWord;
      unsigned p1=0; while(!onOff[p1])++p1;
      unsigned p2=w.length(); while(!onOff[p2-1]) --p2;
      for(unsigned p=p1;p<p2;++p) 
	if(onOff[p]) lWord += w[p];
	else			 lWord += 'N';
      subLWords.push_back(lWord);
      // update the onOff array for the next subLWord
      if(k!=(count-1)) {
	unsigned mpCount=0; //if(onOff[onOff.size()-1])mpCount++;
	unsigned mp=w.length()-1;
	while (!onOff[mp-1] || onOff[mp]) {
	  if(onOff[mp--]) ++mpCount;
	} 
	onOff[mp-1]=false; onOff[mp]=true;
	for(unsigned mp2=mp+1;((mp2<(mp+mpCount+1))&&(mp2<onOff.size()));++mp2)
	  onOff[mp2]=true;
	for(unsigned mp3=(mp+mpCount+1);mp3<onOff.size();++mp3)
	  onOff[mp3]=false;
      }
    }
  }
  return subLWords;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
strIntMap SubWords::operator()(const std::string& w,unsigned m) const {
  strIntMap mWords;
  int k=w.length();
  if(k<2) 
    WriteIntToLogs()(std::string("word length too short in SubWords: ")+w,k);
  else if (m>= k)
    WriteIntToLogs()(std::string("invalid m value in SubWords: "),m);
  else if (!m)
    WriteIntToLogs()(std::string("zero m value in SubWords: "),m);
  else if(!ValidDNASequence()(w) && !ValidRNASequence()(w))
    WriteToLogs()(std::string("word w not a valid nucleic acid sequence in SubWords: ")+w);
  else 
    for(int i=0;i<(k-m+1);++i) {
      std::string mW= w.substr(i,m);
      if(mWords.count(mW)) ++mWords[mW];
      else                 mWords[mW]=1;
    }
  return mWords;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
bool WCSeqMatch::operator()(const std::string& q,const std::string& s) const { 
  if(q.length()!=s.length()) return false;
  for(unsigned i=0;i<q.length();++i)
    switch(q[i]) {
    case 'c':case 'C': if((s[i]!='c')&&(s[i]!='C')) return false; break;
    case 't':case 'T': if((s[i]!='t')&&(s[i]!='T')) return false; break;
    case 'u':case 'U': if((s[i]!='u')&&(s[i]!='U')) return false; break;
    case 'a':case 'A': if((s[i]!='a')&&(s[i]!='A')) return false; break;
    case 'g':case 'G': if((s[i]!='g')&&(s[i]!='G')) return false; break;
      //
    case 'y':case 'Y': if((s[i]!='c')&&(s[i]!='C')&&(s[i]!='t')&&(s[i]!='T')&&
			  (s[i]!='u')&&(s[i]!='U')&&(s[i]!='y')&&(s[i]!='Y')) return false; break;
    case 'r':case 'R': if((s[i]!='a')&&(s[i]!='A')&&(s[i]!='g')&&(s[i]!='G')&&
			  (s[i]!='r')&&(s[i]!='R')) return false; break;
    case 'w':case 'W': if((s[i]!='a')&&(s[i]!='A')&&(s[i]!='t')&&(s[i]!='T')&&
			  (s[i]!='u')&&(s[i]!='U')&&(s[i]!='w')&&(s[i]!='W')) return false; break;
    case 's':case 'S': if((s[i]!='c')&&(s[i]!='C')&&(s[i]!='g')&&(s[i]!='G')&&
			  (s[i]!='s')&&(s[i]!='S')) return false; break;
    case 'k':case 'K': if((s[i]!='g')&&(s[i]!='G')&&(s[i]!='t')&&(s[i]!='T')&&
			  (s[i]!='u')&&(s[i]!='U')&&(s[i]!='k')&&(s[i]!='K')) return false; break;
    case 'm':case 'M': if((s[i]!='c')&&(s[i]!='C')&&(s[i]!='a')&&(s[i]!='A')&&
			  (s[i]!='m')&&(s[i]!='M')) return false; break;
    //
    case 'b':case 'B': if((s[i]!='g')&&(s[i]!='G')&&(s[i]!='t')&&(s[i]!='T')&&
			  (s[i]!='c')&&(s[i]!='C')&&(s[i]!='u')&&(s[i]!='U')&&
			  (s[i]!='b')&&(s[i]!='B')) return false; break;
    case 'd':case 'D': if((s[i]!='g')&&(s[i]!='G')&&(s[i]!='t')&&(s[i]!='T')&&
			  (s[i]!='a')&&(s[i]!='A')&&(s[i]!='u')&&(s[i]!='U')&&
			  (s[i]!='d')&&(s[i]!='D')) return false; break;
    case 'h':case 'H': if((s[i]!='a')&&(s[i]!='A')&&(s[i]!='t')&&(s[i]!='T')&&
			  (s[i]!='c')&&(s[i]!='C')&&(s[i]!='u')&&(s[i]!='U')&&
			  (s[i]!='h')&&(s[i]!='H')) return false; break;
    case 'v':case 'V': if((s[i]!='g')&&(s[i]!='G')&&(s[i]!='a')&&(s[i]!='A')&&
			  (s[i]!='c')&&(s[i]!='C')&&(s[i]!='v')&&(s[i]!='V')) return false; break;
    // 
    case 'n': case 'N': break;
    default: 
      WriteToLog()(std::string("Invalid character in WCSeqMatch"));
      break;
    }
  return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////
iVec WCSeqSearch::operator()(const std::string& q,const std::string& s) const { 
  iVec hits;
  // build the finite state automaton first
  iMat del;
  for(unsigned i=0;i<=q.length();++i) {
    iVec delRow;
    for(unsigned j=0;j<Allbases.length();++j) {
      int k= q.length(); if((i+1)<k) k=i+1;
      bool cont=true;
      while(cont) {
	std::string Pk= q.substr(0,k),Pqa= q.substr(0,i)+Allbases.substr(j,1);
	if(WCSeqMatch()(Pk,Pqa.substr(Pqa.length()-Pk.length(),Pk.length()))) cont=false;
	else --k;
	if(!k) cont=false;
      }
      //			while(!WCSeqMatch()(q.substr(0,k),q.substr(q.length()-k+1,k-1)+Allbases.substr(j,1))&&(k>0)) --k;
      delRow.push_back(k);
    }
    del.push_back(delRow);
  }
  // now run through the subject sequence, saving hits
  unsigned state=0;
  for(unsigned m=0;m<s.length();++m) {
    state= del[state][Allbases.find(s[m])];
    if(state==q.length()) hits.push_back(m-q.length()+1);
  }
  return iVec(hits);
}
//////////////////////////////////////////////////////////////////////////////////////////////
iVec kVector::operator()(const string& s,int k) const {
  if(s.length()<k) return iVec();
  if(k<=0) return iVec(); 
  if(k>maxKvectorK) k=maxKvectorK;
  int nWords=4;
  for(int i=0;i<k;++i) nWords <<= 2;
  //  int nWords=pow(4,k);
  iVec kVec(nWords+1,int(0));
  for(int i=0;i<s.length()-k;++i) {
    int idx= SeqToIdx()(s.substr(i,k));
    if(i>=0) ++kVec[idx];
    else     ++kVec[nWords];
  }
  return kVec;
}
//////////////////////////////////////////////////////////////////////////////////////////////
iMat kPosVector::operator()(const string& s,int k) const {
  if(s.length()<k) return iMat();
  if(k<=0) return iMat(); 
  if(k>maxKvectorK) k=maxKvectorK;
  int nWords=4;//pow(4,k);
  for(int i=0;i<k;++i) nWords <<=2;
  int avgCount= 1+ (s.length()-k+1)/nWords;
  iMat kpVec(nWords+1,iVec());
  for(int i=0;i<nWords+1;++i) kpVec[i].reserve(avgCount);
  for(int i=0;i<s.length()-k;++i) {
    int idx= SeqToIdx()(s.substr(i,k));
    if(i>=0) kpVec[idx].push_back(i);
    else     kpVec[nWords].push_back(i);
  }
  return kpVec;
}
//////////////////////////////////////////////////////////////////////////////////////////////
iVec wordPeriods::operator()(const string& w,bool princ) const {
  if(w.empty()) return iVec();
  iVec p; p.reserve(w.length()-1);
  for(int i=1;i<w.length();++i) {
    bool per=true;
    for(int j=i;(j<w.length())&&per;++j) per= (w[j]==w[j-i]);
    if(per && (!princ || p.empty())) p.push_back(i);
    else if(per && (i%p[0]))         p.push_back(i);
  }
  return p;
}
//////////////////////////////////////////////////////////////////////////////////////////////
iVec TwoWordPeriods::operator()(const strVec& w,bool princ) const {
  if(w.size()<2) return iVec();
  if(w[0].empty()||w[1].empty()) return iVec();
  if(w[0]==w[1]) return wordPeriods()(w[0],princ);
  iVec p; p.reserve(w[0].length()-1);
  for(int i=1;i<w[0].length();++i) {
    bool per=true;
    for(int j=i;(j<w[0].length())&&((j-i)<w[1].length())&&per;++j) per= (w[0][j]==w[1][j-i]);
    if(per && (!princ || p.empty())) p.push_back(i);
    else if(per && (i%p[0]))         p.push_back(i);
  }
  return p;
}
//////////////////////////////////////////////////////////////////////////////////////////////
typedef struct { char base; int n; double p; } baseNpair;
bool operator < (const baseNpair& a,const baseNpair& b) { return (a.p > b.p); }
typedef vector<baseNpair> baseNVec;
//////////////////////////////////////////////////////////////////////////////////////////////
string SeqLogoString::operator()(const dVec& p,int N) const {
  if((p.size()!=4)&&(p.size())!=20) return string();
  bool isProt= (p.size()==20);
  double IC= log2(p.size());
  for(unsigned i=0;i<p.size();++i) IC += p[i]*log2(p[i]);
  int nChars= static_cast<int>(static_cast<double>(N)*IC/log2(static_cast<double>(p.size())));
  baseNVec b;
  for(unsigned i=0;i<p.size();++i) {
    baseNpair temp; 
    if(isProt) temp.base= ProteinBases[i];
    else       temp.base= DNAbases[i];
    temp.n= static_cast<int>(round(nChars*p[i]));
    temp.p= p[i];
    b.push_back(temp);
  }
  std::sort(b.begin(),b.end());
  char buff[501]; std::strstream line(buff,500,std::ios::out);
  if(IC<0.01)     line << "0.00\t";
  else if(IC<0.1) line << std::fixed << std::showpoint << std::setprecision(1) << IC << "\t";
  else            line << std::fixed << std::showpoint << std::setprecision(2) << IC << "\t";
  for(unsigned i=0;i<b.size();++i)
    for(unsigned j=0;j<b[i].n;++j) line << b[i].base;
  line << '\0';
  return buff;
};
//////////////////////////////////////////////////////////////////////////////////////////////
TranslationTable::TranslationTable(const std::string& fname):tMap_(),label_(),desc_() {
  std::ifstream tfile(fname.c_str(),std::ios::in);
  if(tfile.fail()) return;
  // read in the file- expected format is 3 tab-delimited fields, codon, 1-letter, 3-letter
  // also allow option of 2 lines with 2 fields, DESC<tab>description or LABEL<tab>label 
  // if the label or description are not included, the internal data is left blank
  // restrictions on the 3-field lines: 
  //   first field must be three letters, and can only include a,t,u,g, or c; 
  //      case doesn't matter and u's will be converted to t's
  //   second field must be one letter
  //   third field must be three letters
  // we also expect exactly 64 3-field lines with unique first entries
  // errors will be generated if it doesn't follow this format
  char *lineIn= new char[10001];
  int codonCount=0;
  while(!tfile.eof()) {
    tfile.getline(lineIn,10000);
    strVec f= split()(std::string(lineIn),'\t');
    if(f.size()==2) {
      if(f[0]==std::string("DESC")) desc_=f[1];
      else if (f[0]==std::string("LABEL")) label_=f[1];
      else {
	WriteToErr()(std::string("Unknown label in translation file: ")+fname);
	WriteToErr()(std::string(lineIn));
      }
    } else if (f.size()==3) {
      if(f[0].length()!=3 || f[1].length()!=1 || f[2].length()!=3) {
	WriteToErr()(std::string("incorrect field size(s) in translation file: ")+fname);
	WriteToErr()(std::string(lineIn));
      } else {
	string codon= toUpperNoWS()(f[0]);
	if(ValidRNASequence()(codon)) for(int i=0;i<codon.length();++i) if(codon[i]=='U') codon[i]='T';
	if(!ValidDNASequence()(codon)) 
	  WriteToErr()(std::string("illegal codon sequence: ")+codon);
	else if(tMap_.count(codon)) 
	  WriteToErr()(std::string("duplicate codon sequence: ")+codon);
	else {
	  tMap_[codon]= f[1]+f[2]; ++codonCount;
	}
      }
    } else {
      WriteToErr()(std::string("illegal number of fields in translation file: ")+fname);
      WriteToErr()(std::string(lineIn));
    }
  }
  if(codonCount!=64) {
    WriteIntToErr()(std::string("wrong number of entries in translation file: ")+fname+std::string(", "),codonCount);
    tMap_.clear();
  } else 
    WriteToLog()(std::string("successfully loaded translation table from: ")+fname);
  //
  delete[] lineIn;
  tfile.close();
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////
std::string TranslationTable::Translate(const std::string& rna,int frame,bool terse) const {
  if(tMap_.empty()) {
    WriteToErr()(std::string("attempt to use empty Translation table"));
    return std::string();
  }
  if(frame < -3 || frame > 3 || frame==0) {
    WriteIntToErr()(std::string("illegal frame value in TranslationTable::Translate: "),frame);
    return std::string();
  }
  std::string rc;
  if(frame<0) rc= ReverseComplement()(rna);
  const std::string& seq= (frame<0)?rc:rna;
  int start= static_cast<int>(abs(frame)) % 3;
  std::string t;
  int s0=0,s1=1;
  if(!terse) { s0=1; s1=3; }
  while(start+2 < seq.length()) {
    tjCodonMap::const_iterator cdn= tMap_.find(seq.substr(start,3));
    if(cdn!= tMap_.end()) {
      t+= (*cdn).second.substr(s0,s1);
      if(!terse) t+= std::string(" ");
    } else {
      WriteToErr()(std::string("codon not found: ")+seq.substr(start,3));
      if(terse) t += std::string("X");
      else      t += std::string("Unk ");
    }
    start+=3;
  }
  return t;
}
//////////////////////////////////////////////////////////////////////////////////////////////
