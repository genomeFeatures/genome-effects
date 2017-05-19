////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <strstream>
#include <iomanip>
#include "../lib/align.h"
#include "../lib/logFile.h"
#include "../lib/basicSeq.h"
//////////////////////////////////////////////////////////////////////////////////////////////
TJalign::TJalign(const string& s1,const string& s2,int match,int mis,int gInit,int gExt):
  s1_(s1),s2_(s2),match_(match),miss_(mis),gInit_(gInit),gExt_(gExt),S_() {
}
//////////////////////////////////////////////////////////////////////////////////////////////
TJalign::TJalign(const string& s1,const string& s2,const iMat& S,int gInit,int gExt):
  s1_(s1),s2_(s2),match_(),miss_(),gInit_(gInit),gExt_(gExt),S_(S) {
}
//////////////////////////////////////////////////////////////////////////////////////////////
SmithWaterman::SmithWaterman(const string& s1,const string& s2,int match,int mis,int gInit,int gExt):
  TJalign(s1,s2,match,mis,gInit,gExt) {
  S_.erase(S_.begin(),S_.end());
  for(unsigned i=0;i<4;i++) {
    iVec row;
    for(unsigned j=0;j<4;j++)
      if(j==i) row.push_back(match_);
      else     row.push_back(miss_);
    S_.push_back(row);
  }
  BuildMatrixes_();
}
//////////////////////////////////////////////////////////////////////////////////////////////
SmithWaterman::SmithWaterman(const string& s1,const string& s2,const iMat& S,int gInit,int gExt):
  TJalign(s1,s2,S,gInit,gExt) {
  BuildMatrixes_();
}
//////////////////////////////////////////////////////////////////////////////////////////////
int dZero=0, dDown=1, dDiag=2, dAcross=3, dInit=4, dExt=5;
//////////////////////////////////////////////////////////////////////////////////////////////
bool SmithWaterman::BuildMatrixes_(){
  /////////////////////////////////////////////////
  WriteIntToLog()(string("match: "),match_);
  WriteIntToLog()(string("miss : "),miss_);
  WriteIntToLog()(string("init : "),gInit_);
  WriteIntToLog()(string("ext  : "),gExt_);
  WriteToLog()(string("S matrix"));
  for(unsigned i=0;i<S_.size();i++) {
    char line[101]; std::strstream l(line,100); for(unsigned j=0;j<S_[i].size();j++) l << "  " << S_[i][j];
    l << '\0';
    WriteToLog()(string(line));
  }
  WriteToLog()(string("Sequence 1: ")+s1_);
  WriteToLog()(string("Sequence 2: ")+s2_);
  //////////////////////////////////////////////////
  iVec zeroVec; zeroVec.reserve(s2_.length());
  for(unsigned t=0;t<=s2_.length();t++) zeroVec.push_back(0);
  iMat PathL,PathE,PathF,L,E,F;
  PathL.reserve(s1_.length()); PathE.reserve(s1_.length()); PathF.reserve(s1_.length());
  L.reserve(s1_.length()); E.reserve(s1_.length()); F.reserve(s1_.length());
  for(unsigned p=0;p<=s1_.length();p++) {
    PathL.push_back(zeroVec); L.push_back(zeroVec);
    PathE.push_back(zeroVec); E.push_back(zeroVec);
    PathF.push_back(zeroVec); F.push_back(zeroVec);
  }
  // fill the SM array
  //WriteToLog()(string("build score start"));
  int maxP=0,maxT=0; score_=0;
  for(unsigned p=0;p<s1_.length();p++) // current element is always p+1,t+1
    for(unsigned t=0;t<s2_.length();t++) {
      // current element of E (across: insertion in s2, deletion in s1)
      int Einit= L[p+1][t]+gInit_, Eext=  E[p+1][t]+gExt_;
      if((Einit>Eext)&&(Einit>0)) { E[p+1][t+1]= Einit; PathE[p+1][t+1]=dInit; }
      else if(Eext>0)             { E[p+1][t+1]= Eext;  PathE[p+1][t+1]=dExt; }
      else                        PathE[p+1][t+1]=dZero;
      // current element of F (down: insertion in s1, deletion in s2)
      int Finit= L[p][t+1]+gInit_, Fext=  F[p][t+1]+gExt_;
      if((Finit>Fext)&&(Finit>0)) { F[p+1][t+1]= Finit; PathF[p+1][t+1]=dInit; }
      else if(Fext>0)             { F[p+1][t+1]= Fext;  PathF[p+1][t+1]=dExt; }
      else                        PathF[p+1][t+1]=dZero;
      // now L
      int ind1= (S_.size()==4)?DNAbases.find(s1_[p]):ProteinBases.find(s1_[p]);
      int ind2= (S_.size()==4)?DNAbases.find(s2_[t]):ProteinBases.find(s2_[t]);
      int diag=0; bool badIdx=false;
      if((ind1<0)||(((S_.size()==4)&&(ind1>4))||(ind1>23))) {
	WriteIntToErr()(string("out of range base in s1: ")+s1_.substr(p,1)+string(": "),p);
	//WriteToLog()(s1_); 
	badIdx=true;
      }
      if((ind2<0)||(((S_.size()==4)&&(ind2>4))||(ind2>23))) {
	WriteIntToErr()(string("out of range base in s2: ")+s2_.substr(t,1)+string(": "),t);
	//WriteToLog()(s2_); 
	badIdx=true;
      }
      if(!badIdx) diag  = L[p][t]+S_[ind1][ind2];
      else diag=0;
      //      char line[301]; std::strstream buff(line,300,ios::out);
      //      buff << "p: " << p << "; t: " << t << "; ind1: " << ind1 << "; ind2: " << ind2; 
      //      buff << "; diag: " << diag << "; E: " << E[p+1][t+1] << "; F: " << F[p+1][t+1] << '\0';
      //      WriteToLog()(string(line));
      // choose the largest, unless it's less than zero
      if((diag>=E[p+1][t+1])&&(diag>=F[p+1][t+1])&&(diag>0)) {
	PathL[p+1][t+1]= dDiag; L[p+1][t+1]= diag;
      } else if((E[p+1][t+1]>=diag)&&(E[p+1][t+1]>=F[p+1][t+1])&&(E[p+1][t+1]>0)) {
	PathL[p+1][t+1]= dAcross; L[p+1][t+1]= E[p+1][t+1];
      } else if((F[p+1][t+1]>=E[p+1][t+1])&&(F[p+1][t+1]>=diag)&&(F[p+1][t+1]>0)) {
	PathL[p+1][t+1]= dDown; L[p+1][t+1]= F[p+1][t+1];
      } else
	PathL[p+1][t+1]= dZero;
      if(L[p+1][t+1]>score_) {score_=L[p+1][t+1]; maxT=t+1; maxP=p+1;}
    }
  // traceback
  //WriteToLog()(string("matrix built, traceback start"));
  if(score_>0) {
    stop2_=maxT; stop1_=maxP;
    while(L[maxP][maxT]) {
      if(PathL[maxP][maxT]==dDown) { 
	while((PathF[maxP][maxT]==dExt)&&(L[maxP][maxT])) {
	  a1_ += s1_[--maxP]; a2_ += '-';
	}
	a1_ += (s1_[--maxP]+('a'-'A')); a2_ += '-';
      } else if(PathL[maxP][maxT]==dDiag) {
	if(s1_[maxP-1]==s2_[maxT-1]) {
	  a1_ += s1_[--maxP]; a2_ += s2_[--maxT];
	} else {
	  a1_ += (s1_[--maxP]+('a'-'A')); a2_ += (s2_[--maxT]+('a'-'A'));
	}
      } else {
	while((PathE[maxP][maxT]==dExt)&&(L[maxP][maxT])) {
	  a1_ += '-'; a2_ += (s2_[--maxT]+('a'-'A'));
	}
	a1_ += '-'; a2_ += (s2_[--maxT]+('a'-'A'));
      }
    }
    start1_=maxP+1; start2_=maxT+1;
    // reverse the sequences
    for(unsigned i=0;i<a1_.length()/2;i++) {
      int e=a1_.length()-i-1;
      char temp= a1_[i]; a1_[i]=a1_[e]; a1_[e]=temp;
      temp= a2_[i]; a2_[i]=a2_[e]; a2_[e]=temp;
    } 
  }
  //WriteToLog()(string("traceback complete"));
  return (score_>0);
}
//////////////////////////////////////////////////////////////////////////////////////////////
NeedlemanWunsch::NeedlemanWunsch(const string& s1,const string& s2,bool e1,bool e2,int match,int mis,\
				 int gInit,int gExt):TJalign(s1,s2,match,mis,gInit,gExt),pEnd1_(e1),pEnd2_(e2) {
  S_.erase(S_.begin(),S_.end());
  for(unsigned i=0;i<4;i++) {
    iVec row;
    for(unsigned j=0;j<4;j++)
      if(j==i) row.push_back(match_);
      else     row.push_back(miss_);
    S_.push_back(row);
  }
  BuildMatrixes_();
}
//////////////////////////////////////////////////////////////////////////////////////////////
bool NeedlemanWunsch::BuildMatrixes_(){
  iVec zeroVec; zeroVec.reserve(s2_.length());
  for(unsigned t=0;t<=s2_.length();t++) zeroVec.push_back(0);
  iMat PathL,PathE,PathF,L,E,F;
  PathL.reserve(s1_.length()); PathE.reserve(s1_.length()); PathF.reserve(s1_.length());
  L.reserve(s1_.length()); E.reserve(s1_.length()); F.reserve(s1_.length());
  for(unsigned p=0;p<=s1_.length();p++) {
    PathL.push_back(zeroVec); L.push_back(zeroVec);
    PathE.push_back(zeroVec); E.push_back(zeroVec);
    PathF.push_back(zeroVec); F.push_back(zeroVec);
  }
  // fill the SM array
  if(pEnd1_) for(int t=0;t<s2_.length();t++) { L[0][t+1]= gInit_+ t*gExt_; PathL[0][t+1]=dAcross; }
  for(unsigned p=0;p<s1_.length();p++) { // current element is always p+1,t+1
    if(pEnd2_) { L[p+1][0]= gInit_+static_cast<int>(p)*gExt_; PathL[p+1][0]=dDown; }
    for(unsigned t=0;t<s2_.length();t++) {
      // current element of E (across: insertion in s2, deletion in s1)
      if(!pEnd1_ && ((p+1)==s1_.length())) 
	 E[p+1][t+1]=E[p+1][t];
      else {
	int Einit= L[p+1][t]+gInit_, Eext=  E[p+1][t]+gExt_;
	if(Einit>Eext) { E[p+1][t+1]= Einit; PathE[p+1][t+1]=dInit; }
	else           { E[p+1][t+1]= Eext;  PathE[p+1][t+1]=dExt; }
      }
      // current element of F (down: insertion in s1, deletion in s2)
      if(!pEnd2_ && ((t+1)==s2_.length()))
	F[p+1][t+1]=F[p][t+1];
      else {
	int Finit= L[p][t+1]+gInit_, Fext=  F[p][t+1]+gExt_;
	if(Finit>Fext) { F[p+1][t+1]= Finit; PathF[p+1][t+1]=dInit; }
	else           { F[p+1][t+1]= Fext;  PathF[p+1][t+1]=dExt; }
      }
      // now L
      int ind1= (S_.size()==4)?DNAbases.find(s1_[p]):ProteinBases.find(s1_[p]);
      int ind2= (S_.size()==4)?DNAbases.find(s2_[t]):ProteinBases.find(s2_[t]);
      int diag  = L[p][t]+S_[ind1][ind2];
      // choose the largest, unless it's less than zero
      if((diag>=E[p+1][t+1])&&(diag>=F[p+1][t+1])) {
	PathL[p+1][t+1]= dDiag; L[p+1][t+1]= diag;
      } else if((E[p+1][t+1]>=diag)&&(E[p+1][t+1]>=F[p+1][t+1])) {
	PathL[p+1][t+1]= dAcross; L[p+1][t+1]= E[p+1][t+1];
      } else if((F[p+1][t+1]>=E[p+1][t+1])&&(F[p+1][t+1]>=diag)) {
	PathL[p+1][t+1]= dDown; L[p+1][t+1]= F[p+1][t+1];
      } 
    }
  }
  // find the max score, under the constraints of the boundary conditions
  int maxP=s1_.length(),maxT=s2_.length(); 
  if(pEnd2_ && pEnd1_) {
    score_=L[maxP][maxT]; 
  } else if(pEnd2_) {
    score_=-100000;
    for(unsigned j=0;j<L[maxP].size();j++)
      if(L[maxP][j]>score_) { score_=L[maxP][j]; maxT=j; } 
  } else {
    score_=-100000;
    for(unsigned i=0;i<L.size();i++)
      if(L[i][maxT]>score_) { score_=L[i][maxT]; maxP=i; } 
  }
/////////////////////////////////////
  // traceback
  if(!pEnd1_) { // no end penalty on sequence 2
    unsigned tt=s2_.length(); 
    while(tt!=maxT) {
      a1_ += '-'; a2_ += (s2_[--tt]+('a'-'A'));
    }
  }
  if(!pEnd2_) { // no end penalty on sequence 1
    unsigned pp=s1_.length(); 
    while(pp!=maxP) {
      a1_ += (s1_[--pp]+('a'-'A')); a2_ += '-';
    }
  }
  while(maxP||maxT) {
    if(PathL[maxP][maxT]==dDown) { 
      while((PathF[maxP][maxT]==dExt)&&maxP) {
	a1_ += (s1_[--maxP]+('a'-'A')); a2_ += '-';
      }
      if(maxP) { a1_ += (s1_[--maxP]+('a'-'A')); a2_ += '-'; }
    } else if(PathL[maxP][maxT]==dDiag) {
      if(s1_[maxP-1]==s2_[maxT-1]) {
	a1_ += s1_[--maxP]; a2_ += s2_[--maxT];
      } else {
	a1_ += (s1_[--maxP]+('a'-'A')); a2_ += (s2_[--maxT]+('a'-'A'));
      }
    } else {
      while((PathE[maxP][maxT]==dExt)&&maxT) {
	a1_ += '-'; a2_ += (s2_[--maxT]+('a'-'A'));
      }
      if(maxT) { a1_ += '-'; a2_ += (s2_[--maxT]+('a'-'A')); }
    }
  }
  // reverse the sequences
  for(unsigned i=0;i<a1_.length()/2;i++) {
    int e=a1_.length()-i-1;
    char temp= a1_[i]; a1_[i]=a1_[e]; a1_[e]=temp;
    temp= a2_[i]; a2_[i]=a2_[e]; a2_[e]=temp;
  } 
  // if we're not penalizing end gaps, clip them from the sequences
  if(!pEnd1_) {
    int i1=0,i2=a1_.length();
    while(a1_[i1]=='-') ++i1;
    while(a1_[i2-1]=='-') --i2;
    if(i1&&(i2!=a1_.length())) {
      a1_=a1_.substr(i1,i2-i1);
      a2_=a2_.substr(i1,i2-i1);
    }
  }
  if(!pEnd2_) {
    int i1=0,i2=a2_.length();
    while(a2_[i1]=='-') ++i1;
    while(a2_[i2-1]=='-') --i2;
    if(i1&&(i2!=a2_.length())) {
      a1_=a1_.substr(i1,i2-i1);
      a2_=a2_.substr(i1,i2-i1);
    }
  }

  return (score_>0);
}
