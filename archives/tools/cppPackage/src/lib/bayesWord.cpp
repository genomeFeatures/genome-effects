//////////////////////////////////////////////////////////////////////////////////////////////
#include <math.h>
#ifdef INTEL
#include <mathimf.h>
#endif
#include "../lib/bayesWord.h"
#include "../lib/basicSeq.h"
#include "../lib/logFile.h"
#include "../lib/util.h"
#include <algorithm>
static const double errLogVal=-1.0e10;
static const double errLinVal=0.0;
//////////////////////////////////////////////////////////////////////////////////////////////
BayesWordModel::BayesWordModel(const strVec& s,const iVec& p,int N,double pScW,double pScBG,
			       const dVec& totalCounts,int omit,const dMat& prior,SeqType seqType):
 bFr_(totalCounts),bMod_(),wFr_(),wMod_(),bPs_(0.1),wPs_(0.1),seqT_(seqType),N_(N),nHit_(0) {
  if((pScW>=0.0)&&(pScW<=1.0))   wPs_= pScW;
  if((pScBG>=0.0)&&(pScBG<=1.0)) bPs_= pScBG;
  int nBases= (seqT_==Protein)?20:4;
  // other error checks are still needed
  if(s.size()!=p.size()) { 
    WriteToLogs()(string("sequence and position vectors are not same size in BayesWordModel constructor"));
    return;
  }
  if(!bFr_.empty() && bFr_.size()!=nBases) { 
    WriteToLogs()(string("totalCounts vector dimension is wrong in BayesWordModel constructor"));
    return;
  }
  bool withPrior= (prior.size()==N_);
  if(withPrior) {
    for(unsigned i=0;(i<N_)&&withPrior;++i) withPrior= (prior[i].size()==nBases);
  }
  // set up the matrixes first
  // bFr_ is a one dimensional vector of doubles, with nBases elements- initially total count of each basetype
  if(bFr_.empty()) {
    bFr_= dVec(4,double(0.0));
    for(strVec::const_iterator svIt=s.begin();svIt!=s.end();++svIt)
      for(string::const_iterator sIt=(*svIt).begin();sIt!=(*svIt).end();++sIt) { 
	int idx= BaseToIndex(*sIt);
	if(idx>=0) ++bFr_[idx];
      }
  }
  // wFr_ is a two dimensional matrix of doubles, N x nBases, with count of each basetype at each position
  wFr_.reserve(N_); for(unsigned i=0;i<N_;++i) { dVec temp(nBases,double(0.0)); wFr_.push_back(temp);  }
  // process each sequence, building the count matrix, 
  // not including the sequence if either indicated by omit or by p[i]== -1, which means to ignore
  double totalBGCount=0.0,totalModelCount=0.0,totalPriorCount=0.0;
  for(unsigned i=0;i<nBases;++i) totalBGCount+= bFr_[i];
  for(unsigned i=0;i<s.size();++i) 
    if((i!=omit) && (p[i]>=0)) {
      ++nHit_;
      for(unsigned j=p[i];j<(p[i]+N_);++j) 
	if((j>=0)&&(j<s[i].length())) {
	  int idx= BaseToIndex(s[i][j]);
	  if(idx>=0) { ++wFr_[j-p[i]][idx]; ++totalModelCount; --bFr_[idx]; --totalBGCount; }
	}
    }
  if(withPrior) 
    for(unsigned i=0;i<N_;++i) for(unsigned j=0;j<nBases;++j) totalPriorCount+= prior[i][j]; 
  //  WriteFloatToLog()(string("total bg count: "),totalBGCount);
  //  WriteFloatToLog()(string("total model count: "),totalModelCount);
  if(withPrior) WriteFloatToLog()(string("total prior count: "),totalPriorCount);
  // now build the models
  // first background
  bMod_.reserve(nBases);
  double bgInvDenom= 1.0/(totalBGCount+bPs_*static_cast<double>(nBases));
  for(unsigned i=0;i<nBases;++i) bMod_.push_back(log2((bFr_[i]+pScBG)*bgInvDenom));
  // now the word
  double bgRatio= totalModelCount/totalBGCount,priorRatio=(withPrior?totalModelCount/totalPriorCount:0.0);
  double coeff= (withPrior?wPs_*priorRatio:wPs_*bgRatio);
  wMod_.reserve(N_);
  for(unsigned i=0;i<N_;++i) {
    double tPCount=0.0,tPPCount=0.0;
    for(unsigned j=0;j<nBases;++j) { tPCount += wFr_[i][j]; if(withPrior) tPPCount+= prior[i][j]; }
    dVec temp; temp.reserve(nBases);
    double        invDenom= 1.0/(tPCount+coeff*(totalBGCount+static_cast<double>(nBases)));
    if(withPrior) invDenom= 1.0/(tPCount+coeff*(tPPCount+static_cast<double>(nBases)));
    for(unsigned j=0;j<nBases;++j)
      if(withPrior) temp.push_back(log2((wFr_[i][j]+coeff*(prior[i][j]+1.0))*invDenom));
      else          temp.push_back(log2((wFr_[i][j]+coeff*(bFr_[j]+1.0))*invDenom));
    wMod_.push_back(temp); 
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////
BayesWordModel::BayesWordModel(const dMat& wCounts,const dVec& bgCounts,int n,double pScW,double pScBG,
			       SeqType seqType):wFr_(wCounts),bFr_(bgCounts),N_(n),bPs_(0.1),wPs_(0.1),seqT_(seqType) {
  if((pScW>=0.0)&&(pScW<=1.0))   wPs_= pScW;
  if((pScBG>=0.0)&&(pScBG<=1.0)) bPs_= pScBG;
  int nBases= (seqT_==Protein)?20:4;
  if(bFr_.size()!=nBases) return;
  for(int i=0;i<wFr_.size();++i) if(wFr_[i].size()!=nBases) return;
  double TotalWCounts=0.0,TotalBGCounts=0.0;
  for(int i=0;i<nBases;++i) {
      TotalBGCounts += bFr_[i];
      for(int j=0;j<wFr_.size();++j) TotalWCounts+= wFr_[j][i];
  }
  bMod_.reserve(nBases);
  double bgInvDenom= 1.0/(TotalBGCounts+bPs_*static_cast<double>(nBases));
  for(unsigned i=0;i<nBases;++i) bMod_.push_back(log2((bFr_[i]+pScBG)*bgInvDenom));
  double bgRatio= TotalWCounts/TotalBGCounts;
  double coeff= wPs_*bgRatio;
  wMod_.reserve(N_);
  for(unsigned i=0;i<N_;++i) {
    double tPCount=0.0,tPPCount=0.0;
    for(unsigned j=0;j<nBases;++j) tPCount += wFr_[i][j]; 
    dVec temp; temp.reserve(nBases);
    double invDenom= 1.0/(tPCount+coeff*(TotalBGCounts+static_cast<double>(nBases)));
    for(unsigned j=0;j<nBases;++j)
      temp.push_back(log2((wFr_[i][j]+coeff*(bFr_[j]+1.0))*invDenom));
    wMod_.push_back(temp); 
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////
BayesWordModel::BayesWordModel(const dMat& wordModel,const dVec& bgModel,int n,SeqType seqType):
    bFr_(bgModel),bMod_(bgModel),wFr_(wordModel),wMod_(wordModel),bPs_(0.0),wPs_(0.0),seqT_(seqType),N_(n) {
    for(int i=0;i<bFr_.size();++i) bFr_[i]=pow(2.0,bMod_[i]);
    for(int i=0;i<wMod_.size();++i)
	for(int j=0;j<wMod_[i].size();++j)
	    wFr_[i][j]= pow(2.0,wMod_[i][j]);
    return;
}
//////////////////////////////////////////////////////////////////////////////////////////////
inline int BayesWordModel::BaseToIndex(char b) const {
  int idx;
  if     (seqT_==DNA) { idx= DNAbases.find(b);     if(idx>=4) idx= -1; } 
  else if(seqT_==RNA) { idx= RNAbases.find(b);     if(idx>=4) idx= -1; } 
  else                { idx= ProteinBases.find(b); if(idx>=19) idx= -1; }
  return idx;
} 
//////////////////////////////////////////////////////////////////////////////////////////////
double BayesWordModel::logPword(const string& w)    const {
  if(w.length()!=N_) return errLogVal;
  double logp=0.0;
  for(unsigned i=0;i<w.length();++i) {
    int idx= BaseToIndex(w[i]);
    if(idx<0) return errLogVal;
    logp+= wMod_[i][idx];
  }
  return logp;
}
//////////////////////////////////////////////////////////////////////////////////////////////
double BayesWordModel::logPwordRatio(const string& w)    const {
  if(w.length()!=N_) return errLogVal;
  double logp=0.0;
  for(unsigned i=0;i<w.length();++i) {
    int idx=BaseToIndex(w[i]);
    if(idx<0) return errLogVal;
    logp+= (wMod_[i][idx]-bMod_[idx]);
  }
  return logp;
}
//////////////////////////////////////////////////////////////////////////////////////////////
double BayesWordModel::logPwordBG(const string& w)  const {
  if(w.length()!=N_) return errLogVal;
  double logp=0.0;
  for(unsigned i=0;i<w.length();++i) {
    int idx=BaseToIndex(w[i]);
    if(idx<0) return errLogVal;
    logp+= bMod_[idx];
  }
  return logp;
}
//////////////////////////////////////////////////////////////////////////////////////////////
double BayesWordModel::logOddsWord(const string& w) const { return logPwordRatio(w); }
//////////////////////////////////////////////////////////////////////////////////////////////
double BayesWordModel::MAP() const {
  // need total number of pseudocounts and fraction of each type of residue (assume const across model)
  double totalModelCount=0.0,totalBGCount=0.0;
  dVec posModelCount,baseModelCount; int nBases= ((seqT_==Protein)?20:4);
  for(unsigned i=0;i<nBases;++i) { baseModelCount.push_back(0.0); totalBGCount += bFr_[i]; }
  for(unsigned i=0;i<N_;++i)     posModelCount.push_back(0.0);
  for(unsigned i=0;i<N_;++i)
    for(unsigned j=0;j<nBases;++j) {
      baseModelCount[j] += wFr_[i][j];
      posModelCount[i]  += wFr_[i][j];
      totalModelCount   += wFr_[i][j];
    }
  double psRatio= wPs_*totalModelCount/totalBGCount;  
  double totalPC= wPs_*totalModelCount;
  dVec PC; for(unsigned i=0;i<nBases;++i) PC.push_back(psRatio*bFr_[i]);
  // now, compute the MAP value
  // first, the position independent part
  double MAP_= gammln()(totalPC);
  for(unsigned i=0;i<nBases;++i) MAP_ -= gammln()(PC[i]);
  MAP_ *= N_;
  // next, across the positions of the models
  for(unsigned j=0;j<N_;++j) {
    MAP_ -= gammln()(posModelCount[j]+totalPC);
    for(unsigned i=0;i<nBases;++i) MAP_ += gammln()(wFr_[j][i]+PC[i]);
  }
  // finally, normalize to the null model
  for(unsigned i=0;i<nBases;++i)
    MAP_ -= baseModelCount[i]*log(pow(2.0,bMod_[i]));
  return MAP_;
}
//////////////////////////////////////////////////////////////////////////////////////////////
dVec BayesWordModel::SeqProbVector(const string& s) const {
  dVec prob; prob.reserve(s.length());
  if(s.length()>=N_) {
    for(int i=0;i<(s.length()-N_+1);++i) prob.push_back(pow(2.0,logPword(s.substr(i,N_))));
    for(int j=s.length()-N_+1;j<s.length();++j) prob.push_back(0.0);
  }
  return prob;
}
//////////////////////////////////////////////////////////////////////////////////////////////
dVec BayesWordModel::SeqProbRatioVector(const string& s) const {
  dVec prob; prob.reserve(s.length());
  if(s.length()>=N_) {
    double max=0.0;
    for(int i=0;i<(s.length()-N_+1);++i) {
      prob.push_back(pow(2.0,logPwordRatio(s.substr(i,N_))));
      if(prob[i]>max) max=prob[i];
    }
    for(int j=s.length()-N_+1;j<s.length();++j) prob.push_back(0.0);
    for(int i=0;i<prob.size();++i) prob[i] /= max;
  }
  return prob;
}
//////////////////////////////////////////////////////////////////////////////////////////////
dVec BayesWordModel::SeqProbVector(const string& s,double& pMiss) const {
  dVec prob; prob.reserve(s.length()); pMiss=1.0;
  if(!s.empty()) {
    int i=0; double sumPm=0.0;
    for( ;i<s.length()-N_+1;++i) {
      string w(s.substr(i,N_));
      prob.push_back(pow(2.0,logPword(w)));
      sumPm += prob[i];
    }
    while(i++ < s.length()) prob.push_back(0.0);
    pMiss = 1.0/(1.0+2.0*sumPm/static_cast<double>(s.length()-N_));
  }
  return prob;
}
//////////////////////////////////////////////////////////////////////////////////////////////
dVec BayesWordModel::SeqProbRatioVector(const string& s,double& pMiss) const {
  dVec prob; prob.reserve(s.length()); 
  if(!s.empty()) {
    // get the max possible value of logPword
    double maxW=0.0;
    for(int j=0;j<wMod_.size();++j) {
      double maxPos= -1e6;
      for(int k=0;k<wMod_[j].size();++k) if(wMod_[j][k]>maxPos) maxPos=wMod_[j][k];
      maxW+= maxPos;
    }
    maxW= pow(2.0,maxW);
    // now while computing the probability ratio at each position, also find the best match (maxP)
    int i=0; double maxP=-1.0e6;
    for( ;i<s.length()-N_+1;++i) {
      string w(s.substr(i,N_));
      prob.push_back(pow(2.0,logPwordRatio(w)));
      double pW= pow(2.0,logPword(w));
      if(pW>maxP) maxP=pW;
    }
    while(i++ < s.length()) prob.push_back(0.0);
    //    for(int i=0;i<prob.size();++i) prob[i] /= maxP;
    // heuristic- probability of the sequence not including our model, that is, the probability
    // that the whole thing was generated randomly; use the ratio of the probability of the best
    // match in the sequence to the best possible match to our model; the pivot point of 0.05 
    // was developed emprically
    if(2.0*abs(maxP-maxW)/(maxP+maxW)<0.0001) maxP=maxW;
    double ratio=maxP/maxW,pivot=0.1;
    if(ratio>=pivot) pMiss= pivot*(1.0-ratio);
    else          pMiss= 1.0-ratio/pivot;
  }
  return prob;
}
//////////////////////////////////////////////////////////////////////////////////////////////
dVec BayesWordModel::TotalCounts() const {
  dVec tc= bFr_;
  for(unsigned i=0;i<wFr_.size();++i)
    for(unsigned j=0;j<wFr_[i].size();++j)
      tc[j]+=wFr_[i][j];
  return tc;
}
//////////////////////////////////////////////////////////////////////////////////////////////
dVec BayesWordModel::SummedProbVector(const strVec& s) const {
  if(!s.size()) return dVec();
  dVec prob= SeqProbRatioVector(s[0]);
  double maxVal=0.0;
  for(unsigned i=1;i<s.size();++i) {
    dVec prob2= SeqProbRatioVector(s[i]);
    for(unsigned j=0;(j<prob2.size())&&(j<prob.size());++j) {
      prob[j] += prob2[j];
      if(prob[j]>maxVal) maxVal=prob[j];
    }
    for(unsigned j=prob.size();j<prob2.size();++j) {
      prob.push_back(prob2[j]);
      if(prob[j]>maxVal) maxVal=prob[j];
    }
  }
  if(maxVal>0.0) for(unsigned i=0;i<prob.size();++i) prob[i] *= (100.0/maxVal);
  return prob;
}
