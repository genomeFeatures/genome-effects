///////////////////////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <utility>
#include <strstream>
#include "../lib/genome.h"
#include "../lib/direct.h"
#include "../lib/logFile.h"
#include "../lib/basicSeq.h"
///////////////////////////////////////////////////////////////////////////////////////////////////
/*
 protected:
  intRange r_;
  string chr_,org_;
  bool str_;
 public:
  //
  string Sequence() const;
  string RevComp() const;
  //
  bool merge(const genomeSegment& gs);
  bool overlap(const genomeSegment& gs,bool forceStrand=false) const;
  bool contains(const genomeSegment& gs,bool forceStrand=false) const;
 */
///////////////////////////////////////////////////////////////////////////////////////////////////
genomeSegment::genomeSegment(const string& org,const string& chr,int i0,int i1):chr_(chr),org_(org),str_(),r_(make_IntRange(-1,-1)) {
  directory d("/seq/genome/");
  //  if(!d.hasSubDirectory(org)) {
  if(!d.hasEntry(org)) {
    //  strVec sd= d.SubDirectoryList(); for(strVec::iterator f=sd.begin();f!=sd.end();++f) WriteToLog()(string("sd: ")+ *f);
    //  if(std::find(sd.begin(),sd.end(),org_)==sd.end()) {
    WriteToLogs()(string("genomeSegment constructor error: organism not found in /seq/genome/: ")+org_);
    return;
  }
  if(i0<i1) { 
    str_=false; r_= make_IntRange(i0,i1);
  } else {
    str_=true;  r_= make_IntRange(i1,i0);
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
genomeSegment::genomeSegment(const string& org,const string& chr,intRange r,bool str):chr_(chr),org_(org),str_(str),r_(r) {
  directory d("/seq/genome/");
  //  strVec sd= d.SubDirectoryList();
  //  if(std::find(sd.begin(),sd.end(),org_)==sd.end()) {
  //  if(!d.hasSubDirectory(org)) {
  if(!d.hasEntry(org)) {
    r_.first= r_.second= -1;
    WriteToLogs()(string("genomeSegment constructor error: organism not found in /seq/genome/: ")+org_);
    return;
  }
  return;
}    
//////
bool test() {
  return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
string genomeSegment::Sequence() const {
  string dname(string("/seq/genome/")+org_+string("/"));
  directory d(dname);
  string ffname,ifname;
  //  if(d.hasFile(chr_+string(".fa"))) { 
  if(d.hasEntry(chr_+string(".fa"))) { 
    ffname= dname+chr_+string(".fa"); ifname= dname+chr_+string(".fa.idx"); 
    //  } else if(d.hasFile(string("latestgp.fa"))) {
  } else if(d.hasEntry(string("latestgp.fa"))) {
    ffname= dname+string("latestgp.fa"); ifname= dname+string("latestgp.fa.idx"); 
  } else {
    WriteToLogs()(string("genomeSegment::Sequence ERROR: no matching chromosome name, or latestgp.fa for ")+chr_
		  +string(" in directory ")+dname);
    //    return TJfastaSequence();
    return string("ERROR: Sequence Retrieval- chromosome match not found");
  }
  TJfastaSequence temp(ffname,ifname,chr_,r_.first,r_.second-r_.first);
  if(!str_) return temp.Sequence();
  return ReverseComplement()(temp.Sequence());
  //  return TJfastaSequence(temp.Header()+string("; antisense"),ReverseComplement()(temp.Sequence()));
}
///////////////////////////////////////////////////////////////////////////////////////////////////
string genomeSegment::Junction5Sequence(int internalN,int externalN) const {
  string dname(string("/seq/genome/")+org_+string("/"));
  directory d(dname);
  string ffname,ifname;
  //  if(d.hasFile(chr_+string(".fa"))) { 
  if(d.hasEntry(chr_+string(".fa"))) { 
    ffname= dname+chr_+string(".fa"); ifname= dname+chr_+string(".fa.idx"); 
    //  } else if(d.hasFile(string("latestgp.fa"))) {
  } else if(d.hasEntry(string("latestgp.fa"))) {
    ffname= dname+string("latestgp.fa"); ifname= dname+string("latestgp.fa.idx"); 
  } else {
    WriteToLogs()(string("genomeSegment::Sequence ERROR: no matching chromosome name, or latestgp.fa for ")+chr_
		  +string(" in directory ")+dname);
    //    return TJfastaSequence();
    return string("Error in Sequence Retrieval- chromosome match not found");
  }
  int start;
  if(!str_) start= r_.first-externalN;
  else      start= r_.second-internalN;
  TJfastaSequence temp(ffname,ifname,chr_,start,internalN+externalN);
  if(!str_) return temp.Sequence();
  return ReverseComplement()(temp.Sequence());
  //  return TJfastaSequence(temp.Header()+string("; antisense"),ReverseComplement()(temp.Sequence()));
}
///////////////////////////////////////////////////////////////////////////////////////////////////
string genomeSegment::Junction3Sequence(int internalN,int externalN) const {
  string dname(string("/seq/genome/")+org_+string("/"));
  directory d(dname);
  string ffname,ifname;
  //  if(d.hasFile(chr_+string(".fa"))) { 
  if(d.hasEntry(chr_+string(".fa"))) { 
    ffname= dname+chr_+string(".fa"); ifname= dname+chr_+string(".fa.idx"); 
    //  } else if(d.hasFile(string("latestgp.fa"))) {
  } else if(d.hasEntry(string("latestgp.fa"))) {
    ffname= dname+string("latestgp.fa"); ifname= dname+string("latestgp.fa.idx"); 
  } else {
    WriteToLogs()(string("genomeSegment::Sequence ERROR: no matching chromosome name, or latestgp.fa for ")+chr_
		  +string(" in directory ")+dname);
    //    return TJfastaSequence();
    return string("Error in Sequence Retrieval- chromosome match not found");
  }
  int start;
  if(str_) start= r_.first-externalN;
  else     start= r_.second-internalN;
  TJfastaSequence temp(ffname,ifname,chr_,start,internalN+externalN);
  if(!str_) return temp.Sequence();
  return ReverseComplement()(temp.Sequence());
  //  return TJfastaSequence(temp.Header()+string("; antisense"),ReverseComplement()(temp.Sequence()));
}
///////////////////////////////////////////////////////////////////////////////////////////////////
bool genomeSegment::merge(const genomeSegment& gs) {
  if(gs.org_ != org_ || gs.chr_!=chr_ || gs.str_!=str_) return false;
  return ::merge(r_,gs.r_);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
bool genomeSegment::overlap(const genomeSegment& gs,bool forceStrand) const {
  if(gs.org_ != org_ || gs.chr_!=chr_) return false;
  if(forceStrand && gs.str_!=str_) return false;
  return ::overlap(r_,gs.r_);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
bool genomeSegment::contains(const genomeSegment& gs,bool forceStrand) const {
  if(gs.org_ != org_ || gs.chr_!=chr_) return false;
  if(forceStrand && gs.str_!=str_) return false;
  return ::contains(r_,gs.r_);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
/*
class TU {
 protected:
  string org_,chr_,tID_,gID_,pID_;
  int start_,stop_,utr5L__,utr3L_;
  bool str_,complete5pr_,complete3pr_;
  irVec exons_,introns_;
 public:
  typedef enum{ ensembl,sgd,tair } datatype;
  TU(const string& line,datatype dt=ensembl);
 */
///////////////////////////////////////////////////////////////////////////////////////////////////
TU::TU(const string& org,const string& line,datatype dt):org_(org),utr5L_(0),utr3L_(0),cdsL_(0),tL_(0) {
  strVec f=split()(line,'\t');
  if(dt==ensembl) {
    if(f.size()>=18) {
      //      WriteToLog()(string("\n** Transcript Start **"));
      int cdsStart_,cdsStop_;
      chr_=f[6]; if(chr_.substr(0,4)==string("Chr_")) chr_= f[6].substr(4,f[6].length()-4);
      tID_=f[5]; gID_=f[15]; pID_=f[13];
      start_=atoi(f[7].c_str())-1; stop_=atoi(f[8].c_str());
      str_= (atoi(f[9].c_str())<0);
      if(!str_) {
	cdsStart_=atoi(f[10].c_str())-1; cdsStop_=atoi(f[11].c_str());
      } else {
	cdsStart_=atoi(f[11].c_str()); cdsStop_= atoi(f[10].c_str())-1;
      }
      strVec exStruct= split()(f[16],':');
      int nExons= (exStruct.size()+1)/2; exons_.reserve(nExons); introns_.reserve(nExons);
      int sumPieces=0; for(int i=0;i<exStruct.size();++i) sumPieces+=atoi(exStruct[i].c_str());
      if(sumPieces!=(stop_-start_)) 
	WriteToLogs()(tID_+string(" ERROR: sum of exons and introns not equal to stop-start"));

      //      WriteIntToLog()(string("genome start "),start_);
      //      WriteIntToLog()(string("genome stop "),stop_);
      //      WriteIntToLog()(string("cds start "),cdsStart_);
      //      WriteIntToLog()(string("cds stop "),cdsStop_);
      //      WriteToLog()(string("exon-intron structure: ")+f[16]);
      if(!str_) {
	int exStart=start_; 
	for(int i=0;i<nExons;++i) {
	  int nextEx= atoi(exStruct[2*i].c_str());
	  exons_.push_back(make_IntRange(exStart,exStart+nextEx)); 
	  //
	  if(cdsStart_>exons_[i]) utr5L_+= ::length(exons_[i]);
	  else if(contains(exons_[i],cdsStart_)) utr5L_+= (cdsStart_-exons_[i].first);
	  if(cdsStop_<exons_[i]) utr3L_+= ::length(exons_[i]);
	  else if(contains(exons_[i],cdsStop_)) utr3L_+= (exons_[i].second-cdsStop_);
	  tL_ += ::length(exons_[i]);
	  //
	  //	  char buff[101]; std::strstream os(buff,100,std::ios::out);  os << "s   exon: " << exStart << "->" << (exStart+nextEx) << '\0';
	  //	  WriteToLog()(string(buff));
	  if((i+1)!=nExons) {
	    int nextInt=atoi(exStruct[2*i+1].c_str());
	    introns_.push_back(make_IntRange(exStart+nextEx,exStart+nextEx+nextInt));
	    //	    os.freeze(false); os.seekp(0,std::ios::beg); 
	    //	    os << "s intron: " << (exStart+nextEx) << "->" << (exStart+nextEx+nextInt) << '\0';
	    //	    WriteToLog()(string(buff));
	    exStart+= nextEx+nextInt;
	  }
	}
      } else {
	int exStop=stop_;
	for(int i=nExons-1;i>=0;--i) {
	  int nextEx= atoi(exStruct[2*i].c_str());
	  exons_.push_back(make_IntRange(exStop-nextEx,exStop));
	  //
	  intRange r= *(exons_.rbegin());
	  if(cdsStart_<r) utr5L_+= ::length(r);
	  else if(contains(r,cdsStart_)) utr5L_+= (r.second-cdsStart_);
	  if(cdsStop_>r) utr3L_+= ::length(r);
	  else if(contains(r,cdsStop_)) utr3L_+= (cdsStop_-r.first);
	  tL_ += ::length(r);
	  //	  
	  //	  WriteIntToLog()(string("temp utr5L "),utr5L_);
	  //	  WriteIntToLog()(string("temp utr3L "),utr3L_);
	  //	  char buff[101]; std::strstream os(buff,100,std::ios::out);  os << "as   exon: " << (exStop-nextEx) << "->" << exStop << '\0';
	  //	  WriteToLog()(string(buff));
	  if(i>0) {
	    int nextInt= atoi(exStruct[2*i-1].c_str());
	    introns_.push_back(make_IntRange(exStop-nextEx-nextInt,exStop-nextEx));
	    //os.freeze(false); os.seekp(0,std::ios::beg); 
	    //os << "as intron: " << (exStop-nextEx-nextInt) << "->" << (exStop-nextEx) << '\0';
	    //WriteToLog()(string(buff));
	    exStop-= (nextEx+nextInt);
	  }
	}
      }
      cdsL_= tL_-utr3L_-utr5L_;
      //      WriteIntToLog()(string("utr5 length "),utr5L_);
      //      WriteIntToLog()(string("utr3 length "),utr3L_);      
      //      WriteToLog()(string("** Transcript Stop **"));
    }
  } else if(dt==sgd) {
  } else if(dt==tair) {
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////
string TU::Sequence() const {
  string s;
  for(int i=0;i<exons_.size();++i) {
    genomeSegment gs(org_,chr_,exons_[i],str_);
    if(gs.begin()<0 || gs.end()<0) return string("ERROR: TU::Sequence, segment has invalid boundaries");
    string es= gs.Sequence();
    if(es.substr(0,5)==string("ERROR")) return string("ERROR: TU::Sequence, segment has invalid sequence name");
    s += es;
  }
  return s;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
gsVec TU::Introns() const {
  gsVec rIntrons; rIntrons.reserve(introns_.size());
  for(int i=0;i<introns_.size();++i) {
    genomeSegment gs(org_,chr_,introns_[i],str_);
    rIntrons.push_back(gs);
  }
  return rIntrons;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
gsVec TU::Exons() const {
  gsVec rExons; rExons.reserve(exons_.size());
  for(int i=0;i<exons_.size();++i) {
    genomeSegment gs(org_,chr_,exons_[i],str_);
    rExons.push_back(gs);
  }
  return rExons;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
strVec TU::IntronSequences() const {
  strVec isv;
  for(int i=0;i<introns_.size();++i) {
    genomeSegment gs(org_,chr_,introns_[i],str_);
    if(gs.begin()<0 || gs.end()<0) return strVec(1,string("ERROR: TU::IntronSequences, segment has invalid boundaries"));
    string is= gs.Sequence();
    if(is.substr(0,5)==string("ERROR")) return strVec(1,string("ERROR: TU::IntronSequences, segment has invalid sequence name"));
    isv.push_back(is);
  }
  return isv;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
strVec TU::ExonSequences() const {
  strVec esv;
  for(int i=0;i<exons_.size();++i) {
    genomeSegment gs(org_,chr_,exons_[i],str_);
    if(gs.begin()<0 || gs.end()<0) return strVec(1,string("ERROR: TU::ExonSequences, segment has invalid boundaries"));
    string es= gs.Sequence();
    if(es.substr(0,5)==string("ERROR")) return strVec(1,string("ERROR: TU::ExonSequences, segment has invalid sequence name"));
    esv.push_back(es);
  }
  return esv;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
strVec TU::ExonIntronSequences(int eN,int iN) const {
  strVec eisv;
  for(int i=0;i<exons_.size()-1;++i) {
    genomeSegment gs(org_,chr_,exons_[i],str_);
    if(gs.begin()<0 || gs.end()<0) return strVec(1,string("ERROR: TU::ExonIntronSequences, segment has invalid boundaries"));
    string es= gs.Junction3Sequence(eN,iN);
    if(es.substr(0,5)==string("ERROR")) return strVec(1,string("ERROR: TU::ExonIntronSequences, segment has invalid sequence name"));
    eisv.push_back(es);
  }
  return eisv;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
strVec TU::IntronExonSequences(int iN,int eN) const {
  strVec iesv;
  for(int i=1;i<exons_.size();++i) {
    genomeSegment gs(org_,chr_,exons_[i],str_);
    if(gs.begin()<0 || gs.end()<0) return strVec(1,string("ERROR: TU::IntronExonSequences, segment has invalid boundaries"));
    string es= gs.Junction5Sequence(eN,iN);
    if(es.substr(0,5)==string("ERROR")) return strVec(1,string("ERROR: TU::IntronExonSequences, segment has invalid sequence name"));
    iesv.push_back(es);
  }
  return iesv;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
int TU::length() const {
  int l=0;
  for(int i=0;i<exons_.size();++i) l+= (exons_[i].second-exons_[i].first);
  return l;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
string TU::UTR5() const {
  string s;
  int utrLength=0; 
  for(int i=0;i<exons_.size()&&s.length()<utr5L_;++i) {
    genomeSegment gs(org_,chr_,exons_[i],str_);
    if(gs.begin()<0 || gs.end()<0) return string("ERROR: TU::Sequence, segment has invalid boundaries");
    string es= gs.Sequence();
    if(es.substr(0,5)==string("ERROR")) return string("ERROR: TU::Sequence, segment has invalid sequence name");
    s += es;
  }
  return s.substr(0,utr5L_);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
string TU::UTR3() const {
  string s;
  for(int i=exons_.size()-1;i>=0 && s.length()<utr3L_;--i) {
    genomeSegment gs(org_,chr_,exons_[i],str_);
    if(gs.begin()<0 || gs.end()<0) return string("ERROR: TU::Sequence, segment has invalid boundaries");
    string es= gs.Sequence();
    if(es.substr(0,5)==string("ERROR")) return string("ERROR: TU::Sequence, segment has invalid sequence name");
    s.insert(0,es);
  }
  return s.substr(s.length()-utr3L_,utr3L_);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
string TU::Upstream(int N) const {
  if(!str_) {
    if(N>start_) N=start_;
    genomeSegment gs(org_,chr_,make_IntRange(start_-N,start_),str_);
    return gs.Sequence();
  } 
  genomeSegment gs(org_,chr_,make_IntRange(stop_,stop_+N),str_);
  return gs.Sequence();
}
///////////////////////////////////////////////////////////////////////////////////////////////////
string TU::Downstream(int N) const {
  if(str_) {
    if(N>start_) N=start_;
    genomeSegment gs(org_,chr_,make_IntRange(start_-N,start_),str_);
    return gs.Sequence();
  }
  genomeSegment gs(org_,chr_,make_IntRange(stop_,stop_+N),str_);
}
/*
  intRange r_;
  string chr_,org_;
  bool str_;

  string rID_,stID_,tID_,detailed_,method_;
 */
///////////////////////////////////////////////////////////////////////////////////////////////////
Repeat::Repeat(const string& line) {
  strVec fields= split()(line);
  if(fields.size()<8) {
    WriteToLogs()(string("ERROR in Repeat(line) constructor, fields<8:  ")+line);
    return;
  } else if(fields.size()>9) 
    WriteToLogs()(string("WARNING in Repeat(line)constructor, fields>9:  ")+line);
  rID_=fields[0]; stID_=fields[1]; tID_=fields[2], chr_=fields[3]; 
  r_= make_IntRange(atoi(fields[4].c_str()),atoi(fields[5].c_str()));
  str_= (*(fields[6].begin())) =='-';
  detailed_= fields[7];
  if(fields.size()>8) method_=fields[8];
}
///////////////////////////////////////////////////////////////////////////////////////////////////
Repeat::Repeat(std::istream& is) {
  char lineIn[1001];
  is.getline(lineIn,1000);
  string line(lineIn);
  strVec fields= split()(line);
  if(fields.size()<8) {
    WriteToLogs()(string("ERROR in Repeat(line) constructor, fields<8:  ")+line);
    return;
  } else if(fields.size()>9) 
    WriteToLogs()(string("WARNING in Repeat(line)constructor, fields>9:  ")+line);
  rID_=fields[0]; stID_=fields[1]; tID_=fields[2], chr_=fields[3]; 
  r_= make_IntRange(atoi(fields[4].c_str()),atoi(fields[5].c_str()));
  str_= (*(fields[6].begin())) =='-';
  detailed_= fields[7];
  if(fields.size()>8) method_=fields[8];
}
///////////////////////////////////////////////////////////////////////////////////////////////////

    
