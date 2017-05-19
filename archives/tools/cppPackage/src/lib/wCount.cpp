#include <algorithm>
#include <math.h>
#ifdef INTEL
#include <mathimf.h>
#endif
#include <time.h>
#include <iostream>
#include <strstream>
#include <fstream>
#include "../lib/logFile.h"
#include "../lib/basicSeq.h"
#include "../lib/wCount.h"
#include "../lib/wordStat.h"
#include "../lib/util.h"
#include "../lib/vectorF.h"
///////////////////////////////////////////////////////////////////////////////////////////
// class cjWCount- a word counting class
// variables:
//  unsigned mN_,nS_; // mN_ = max word length; nS_= number of seqs included in the count
//  uVec nW_;         // (vector of unsigned) nW_[i] = number of words of length i+1
//  siMapVec wc_;     // the counts, a vector of maps from string to integer 
//                        different maps for each size word
//  bool ds_;         // ds_ = true if we are counting double stranded
//
///////////////////////////////////////////////////////////////////////////////////////////
// cjWCount: base constructor, empty counter
cjWCount::cjWCount(unsigned maxN,bool ds):mN_(maxN),nS_(0),nW_(),wc_(),ds_(ds){
  // check that max word size is allowed
  if(mN_>maxWordSize_) {
    WriteIntToLog()(string("max word size too large, resetting: "),maxN);
    mN_=maxWordSize_;
  }
  // initiate count holders
  for(unsigned i=0;i<mN_;i++) {
    strIntMap dummy;
    nW_.push_back(0); wc_.push_back(dummy);
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
// cjWCount: constructor with input sequences
cjWCount::cjWCount(const strVec& seqs,unsigned maxN,bool ds):
mN_(maxN),nS_(0),nW_(),wc_(),ds_(ds) {
  // check max word size
  if(mN_>maxWordSize_) {
    WriteIntToLog()(string("max word size too large, resetting: "),maxN);
    WriteIntToErr()(string("max word size too large, resetting: "),maxN);
    mN_=maxWordSize_;
  }
  // initiate count holders
  for(unsigned i=0;i<mN_;i++) {
    strIntMap dummy;
    nW_.push_back(0); wc_.push_back(dummy);
  }
  // add the input sequences
  AddSequences(seqs);
}
//////////////////////////////////////////////////////////////////////////////////////////////
// assignment operator
cjWCount& cjWCount::operator=(const cjWCount& o){
  if(o== *this) return *this;
  mN_=o.mN_; nS_=o.nS_; nW_=o.nW_; wc_=o.wc_; ds_=o.ds_;
  return *this;
}
///////////////////////////////////////////////////////////////////
// cjWCount Method: AddSequences 
//    adds a vector of sequences; simple loop around calls to AddSequence
void cjWCount::AddSequences(const strVec& seqs){
  for(unsigned s=0;s<seqs.size();s++) 
    AddSequence(seqs[s]);
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount Method: AddSequence
//    adds a sequence to the count matrix. if the sequence is N long, then N mononucleotides
//    are added, N-1 dinucleotides, and so on up to mN_ (max size);  if counting is double 
//    stranded, then words are only stored as one of the complementary pair
void cjWCount::AddSequence(const string& seq){
  nS_++; string rc;
  if(ds_) rc= ReverseComplement()(seq); // faster to do one rc at the beginning 
  double increment=(ds_?2:1);           // for ds counting, add two at a time
  for(unsigned p=0;p<seq.length();p++) 
    for(unsigned w=0;w<mN_;w++) 
      if((p+w)<seq.length()) {
	string word(seq.substr(p,w+1));
	if(ds_) {
	  string rcW=rc.substr(rc.length()-w-1-p,w+1);
	  if(rcW<word) word=rcW;
	}
	if(ValidDNASequence()(word)) {
	  if(wc_[w].count(word)) wc_[w][word]+=increment;
	  else                   wc_[w][word]=increment;
	  nW_[w]+=increment;
	} else if(!wc_[w].count(string("ambig"))) 
	  wc_[w][string("ambig")]=increment;
	else
	  wc_[w][string("ambig")]+=increment;
      }
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount Method: WordCount
//    returns the count for the entered word
unsigned cjWCount::WordCount(const string& w) const {
  unsigned wL=w.length()-1;
  // words longer than mN_ aren't counted, so return 0
  if(wL>=mN_) return 0; 
  // now, check for wildcards
  strVec expWords=ExpandWildCards()(w,'D');
  if(expWords.size()>1) {
    string msg("Expansion of w in cjWCount::WordCount: "); msg+=w;
    msg += string (" -> ");
    for(unsigned j=0;j<expWords.size();j++) 
      if(j) msg += string(",")+expWords[j];
      else  msg += expWords[j];
    WriteToLog()(msg);
  }
  // loop through the expanded words, summing all of the counts
  unsigned count=0;
  for(unsigned i=0;i<expWords.size();i++) { 
    strIntMap::const_iterator iter= wc_[wL].find(expWords[i]);
    if(iter!=wc_[wL].end()) 
      count += (*iter).second;
    else {
      iter= wc_[wL].find(ReverseComplement()(expWords[i]));
      if(iter!=wc_[wL].end()) count += (*iter).second;
    }
  }
  return count;
}
/////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount Method: WordFrac
//    returns the word's fraction of all counted of the same length
double cjWCount::WordFrac(const string& word) const {
  unsigned wL=word.length()-1;
  if(wL>=mN_) return 0;
  unsigned count=WordCount(word);
  if(count) 
    return( static_cast<double>(count)/static_cast<double>(nW_[wL]) );
  else 
    return 0.0;
}
/////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount Method: WordPerc
//    returns the word's percentage of all counted of the same length
double cjWCount::WordPerc(const string& word) const {
  return WordFrac(word)*100.0;
}
/////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount Method: ZScore
//    calculates Z, based on a normal approximation.  As currently implemented, it's only 
//    correct for order = wordlength-2 (4th order for a hexamer)
double cjWCount::ZScore(const string& word,unsigned order) const {
  if(!WordCount(word)) return 0.0;
  if(order>(word.length()-2)) {
    WriteIntToErr()(string("Markov order too high for word: ")+word+string(": "),order);
    WriteIntToLog()(string("Markov order too high for word: ")+word+string(": "),order);
    return -1.0e10;
  }
  //  double E= Expect(word,order);
  unsigned orderP1=order+1;
  double Counts=WordCount(word);
  double num=1.0,den=1.0; dVec hWCounts,gWCounts;
  for(unsigned i=0;i<(word.length()-order);i++) {
    double nc= static_cast<double>(WordCount(word.substr(i,orderP1)));
    num *= nc; hWCounts.push_back(nc);
    if(i) {
      double dc= static_cast<double>(WordCount(word.substr(i,order)));
      den *= dc; gWCounts.push_back(dc);
    }
  }
  // expected value E
  double E;
  if(den!=0.0) E=num/den;
  else         E=0.0;
  // estimated variance V
  double V=E;
  if(den!=0.0) {
    V /= (den*den);
    for(unsigned k=0;k<gWCounts.size();k++) 
      V *= (gWCounts[k]-hWCounts[k])*(gWCounts[k]-hWCounts[k+1]);
  }
//
  double Z=0.0;
  if(V!=0.0) Z= (Counts-E)/sqrt(V);
  return double(Z);
}
////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount Method: KarlinRho
//    calculates the rho score as defined by Karlin and Burge; a cascading ratio that attempts
//    to correct longer words for more significant shorter words
//  06 Aug 2002- this should be re-written to include pseudo counts
double cjWCount::KarlinRho(const string& w) const {
  if(!WordCount(w)||(w.length()<1)) return 1.0;
  if(w.length()==1) return WordPerc(w)/0.25;
  double wp=WordFrac(w);
  for(unsigned l=1;l<w.length();l++) {
    strVec subWords=GetSubLWords()(w,w.length()-l);
    for(unsigned i=0;i<subWords.size();i++) {
      string sw(subWords[i]);
      double ff=WordFrac(sw);
      if(l%2) wp /= ff;//WordFrac(sw);
      else    wp *= ff;//WordFrac(sw);
    }
  }
  return wp;
}
/////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount Method: Entropy
//   calculates Shannon's entropy for the measured counts
//    06 Aug 2002: currently only valid for nucleic acids.  needs to be fixed for proteins
double cjWCount::Entropy(unsigned ord) const {
  unsigned wL=ord;
  if(wL>=mN_) return 1.0e50; // check the order first
  double ent=0.0;
  // we sum over all of the n-mers specified by ord
  unsigned terms=pow(4,ord);
  for(unsigned i=0;i<terms;i++) {
    string source= IdxToSeq()(i,ord); // source is the ith word of length ord
    double fSource=WordFrac(source);  // frequence of the source
    for(unsigned j=0;j<4;j++) {
      string word= source+IdxToSeq()(j,1); // extension- the source plus one base (each possible)
      if(ord) { 
	double wCount=WordCount(word),sCount=WordCount(source);
	// the inclusion of fSource as weight makes the maximum sum log2(alph size), even for 
	// for higher order than 1
	if(wCount>0) ent -= fSource*wCount/sCount*log2(wCount/sCount);
      } else {
	double wFrac=WordFrac(word);
	if(wFrac>0.0) ent -= wFrac*log2(wFrac);
      }
    }
  }
  return ent;
}
////////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount method: Information
//    calculates Shannon's information- 
//    06 Aug 2002: currently only valid for nucleic acids.  needs to be fixed for proteins
double cjWCount::Information(unsigned ord) const {
  unsigned wL=ord-1;
  if(wL>=mN_) return 1.0e50;
  return (2.0-Entropy(ord));
}
////////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount method: AmbiguousCount
//    returns the number of ambiguous counts (words with wildcards, etc) of length N 
unsigned cjWCount::AmbiguousCount(unsigned N) const {
  if(N<=mN_) {
    strIntMap::const_iterator f= wc_[N-1].find(string("ambig"));
    if(f!=wc_[N-1].end()) return (*f).second;
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount method: UniqueWordCount
//    returns the number of different words of length N counted
unsigned cjWCount::UniqueWordCount(unsigned N) const {
  if(N<=mN_) {
    if(wc_[N-1].count(string("ambig"))) return (wc_[N-1].size()-1);
    return wc_[N-1].size();
  } else      return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////
// cJWCount method: TotalWordCount
//    returns the total number of valid words of length N counted
unsigned cjWCount::TotalWordCount(unsigned N) const {
  if(N<=mN_) return nW_[N-1];
  else       return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount Method: Expect
//   returns the expected counts of the word w, based on a o-order Markov prediction
double cjWCount::Expect(const string& w,unsigned o) const {
  if(o>w.length()-2) {
    WriteIntToErr()(string("ExpectFreq: order is too high for word: ")+w,o);
    return double(0.0);
  }
  if(o>(mN_-1)) {
    WriteIntToErr()(string("ExpectFreq: order is too high for maxN: "),o);
    return double(0.0);
  }
  // following section removed 07 Aug 2002- doesn't matter if the word's not there!
  //  strIntMap::const_iterator iter= wc_[w.length()-1].find(w);
  //  if(iter==wc_[w.length()-1].end()) iter= wc_[w.length()-1].find(ReverseComplement()(w));
  //  if(iter==wc_[w.length()-1].end()) return double(0.0);
  //  string word= (*iter).first;
  unsigned orderP1=o+1; //dVec hWCounts,gWCounts; 
  double num=1.0,den=1.0;
  for(unsigned i=0;i<(w.length()-o);i++) {
    num *= static_cast<double>(WordCount(w.substr(i,orderP1)));
    if(i) den *= static_cast<double>(WordCount(w.substr(i,o)));
    //    string hWord(word.substr(i,orderP1));
    //    hWCounts.push_back(static_cast<double>(WordCount(hWord)));
    //    if(i) {
    //      string gWord(word.substr(i,o));
    //      gWCounts.push_back(static_cast<double>(WordCount(gWord)));
    //    }
  }
//  for(unsigned j=0;j<hWCounts.size();j++) num *= hWCounts[j];
//  for(unsigned k=0;k<gWCounts.size();k++) den *= gWCounts[k];
// expected value E
  double E=0.0;
  if(den!=0.0) E=num/den;
//	WriteFloatToLog()(string("exp freq: ")+w,E);
  return E;
}
////////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount Method: ExpectM
//    returns expected number of counts, based on Markov chain of order o, scaled to a new length L
double cjWCount::ExpectM(const string& w,unsigned o,unsigned L) const {
  return (ExpectFreq(w,o)*static_cast<double>(L));
}
////////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount Method: ExpectFreq
//    returns expected number of counts, based on Markov chain of order o, scaled to a new length L
double cjWCount::ExpectFreq(const string& w,unsigned o) const {
// E is the expected count;  to get frequency, normalize to the number of words of the same length
// that would have been counted, but we can't do that directly, since the counts only have to be up
// to size o+1 in order to successfully.  But we know that for every sequence, the difference in the
// counted number of different sizes is just the difference in their lengths-
  return Expect(w,o)/static_cast<double>(nW_[o]-nS_*(w.length()-o-1));
}
////////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount Method: GenerateSequence
//    generates a random sequence, of length N, assuming a Markov chain of the order given
static long seed=0; static unsigned rCount=0;
string cjWCount::GenerateSequence(unsigned N,unsigned order) const {
  if(order>mN_) {
    WriteIntToErr()(string("GenSeq order too high: "),order); return false;
  }
  string gs; gs.reserve(N);
  // generate a new seed at the start, and after 50,000 steps
  if(!seed||(rCount>50000)) { time_t t1; time(&t1); seed= -1*t1; rCount=0; }
  rCount += N;
  while (gs.length()<N) {
    double r= rand3()(&seed);
    double sc=1.0; string source;
    if(order) {
      if(gs.length()>order) source=gs.substr(gs.length()-order,order);
      else                  source=gs;
      if(source.length())   sc=WordFrac(source);
    }
    if(!ds_ || gs.length()) {
      if(r < WordFrac(source+string("C"))/sc)	                                  gs += 'C';
      else if(r < (WordFrac(source+string("C"))+WordFrac(source+string("T")))/sc) gs += 'T';
      else if(r < (WordFrac(source+string("C"))+WordFrac(source+string("T"))+
		   WordFrac(source+string("A")))/sc)			          gs += 'A';
      else                                                                        gs += 'G';
    } else {
      double r2=rand3()(&seed);
      if( (r < WordFrac(string("C")))&&(r2<0.5) ) gs += 'C';
      else if(r < WordFrac(string("C")))          gs += 'G';
      else if(r2<0.5)                             gs += 'A';
      else                                        gs += 'T';
    }
  }
//
  return string(gs);
}
//////////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount Method: ListWordCounts
//    produces a string listing counts- one word and its count per line, tab delimited
string cjWCount::ListWordCounts() const {
  string str;
  for(unsigned i=0;i<mN_;i++) 
    for(strIntMap::const_iterator p=wc_[i].begin();p!=wc_[i].end();p++) { 
      unsigned count= WordCount((*p).first);
      char num[10]; std::strstream numStr(num,10,std::ios::in); numStr << count; //itoa(count,num,10); 
      str += (*p).first;
      if(!ds_) str += (string("/")+ReverseComplement()((*p).first));
      str += (string("\t")+string(num)+string("\n"));
    }
  return string(str);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount Method: ListWords
//    generates a string vector of all words counted.
strVec cjWCount::ListWords() const {
  strVec words;
  for(unsigned i=0;i<mN_;i++) 
    for(strIntMap::const_iterator p=wc_[i].begin();p!=wc_[i].end();p++) 
      if((*p).first!=string("ambig")) {
	words.push_back((*p).first);
      }
  return strVec(words);
}
////////////////////////////////////////////////////////////////////////////////////////////////
// cjWCount Method LoadCountsFromFile
//    reads in counts, assuming ascii tab-delimited word and count, one pair per line
unsigned cjWCount::LoadCountsFromFile(const string& fn) {
  std::ifstream cntFile;
  //  cntFile.open(fn.c_str(),std::ios::in|ios::nocreate); ios::nocreate is deprecated
  cntFile.open(fn.c_str(),std::ios::in);
  if(!(cntFile.fail()||cntFile.eof()) && cntFile) {
    mN_=nS_=0;
    wc_.erase(wc_.begin(),wc_.end());
    nW_.erase(nW_.begin(),nW_.end());
    bool cont=true; unsigned twc=0;
    while(cont&&!cntFile.eof()) {
      string ww; 
      char ch; cntFile.get(ch); 
      while((ch!=' ')&&(ch!='\n')&&(ch!='\t')&&!cntFile.eof()) {
	ww += ch;
	cntFile.get(ch);
      }
      if((ww==string("Ambig-0"))||cntFile.eof()) 
	cont=false;
      else if(ValidDNASequence()(ww)) {
	unsigned cnt; char dummy[1000];
	twc++;
	cntFile >> cnt; cntFile.getline(dummy,999);
	if(ww.length()>mN_) {
	  mN_ = ww.length();
	  strIntMap temp; temp[ww]= cnt;
	  wc_.push_back(temp);
	  nW_.push_back(1);
	} else {
	  wc_[ww.length()-1][ww]= cnt;
	  nW_[ww.length()-1]++;
	}
      } else {
	WriteToErr()(string("invalid word in bg file: ")+ww);
	WriteToLog()(string("invalid word in bg file: ")+ww);
	cont=false;
      }
    }
    return twc;
  }
  return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
// class cjPWCount- derived from cjWCount
//   word counting class that also maintains position information
//    new variables:
//  sivMapVec wp_;              vector of maps from string to vector of ints
//  unsigned minSeqL_,maxSeqL_; unsigned
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
// cjPWCount base constructor
//    constructor to cjWCount will check values of maxN
cjPWCount::cjPWCount(unsigned maxN,bool ds):cjWCount(maxN,ds),wp_(),minSeqL_(0),maxSeqL_(0){
  for(unsigned i=0;i<mN_;i++) {
    strIVecMap dummy;    // empty string to integer vector map
    wp_.push_back(dummy);
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////
// cjPWCount constructor with input data
//   
cjPWCount::cjPWCount(const strVec& seqs,unsigned maxN,bool ds):
  cjWCount(seqs,maxN,ds),wp_(),minSeqL_(0),maxSeqL_(0) {
  for(unsigned i=0;i<mN_;i++) {
    strIVecMap dummy;     // empty string to integer vector map
    wp_.push_back(dummy);
  }
  AddSequences(seqs); // add the sequences
}
/////////////////////////////////////////////////////////////////////////////////////////////////
// cjPWCount Method: AddSequence
//    adds the sequence to the counts
void cjPWCount::AddSequence(const string& seq){
  nS_++;  // increment the number of sequences
  string rc; if(ds_) rc=ReverseComplement()(seq);
  double increment=ds_?2:1;
  if((seq.length()<minSeqL_)||(nS_==1)) minSeqL_=seq.length();
  if((seq.length()>maxSeqL_)||(nS_==1)) maxSeqL_=seq.length();
  for(unsigned p=0;p<seq.length();p++)
    for(unsigned w=0;w<mN_;w++) 
      if((p+w)<seq.length()) {
	string word(seq.substr(p,w+1));
	if(ds_) {
	  string rcW=rc.substr(rc.length()-w-1-p,w+1);
	  if(rcW<word) word=rcW;
	}
	if(ValidDNASequence()(word)) {
	  if(wc_[w].count(word)) {
	    wc_[w][word]+=increment; wp_[w][word].push_back(p);
	  } else {
	    wc_[w][word]=increment; 
	    iVec temp; temp.push_back(p); wp_[w][word]= temp; 
	  }
	  nW_[w]+=increment;
	} else if(!wc_[w].count(string("ambig"))) 
	  wc_[w][string("ambig")]=increment;
	else
	  wc_[w][string("ambig")]+=increment;
	//	WriteToLog()(string("invalid word not counted: ")+word);
      }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////
// cjPWCount assignment operator
cjPWCount& cjPWCount::operator=(const cjPWCount& o){
  if(o== *this) return *this;
  mN_=o.mN_; nS_=o.nS_; nW_=o.nW_; wc_=o.wc_; ds_=o.ds_; 
  wp_=o.wp_; minSeqL_=o.minSeqL_; maxSeqL_=o.maxSeqL_;
  return *this;
}
///////////////////////////////////////////////////////////////////////////////////////////////
// cjPWCount Method: PosMean
//    returns the mean (average) position of the word, averaged over all occurrences
double cjPWCount::PosMean(const string& w) const {
  strIVecMap::const_iterator iter= wp_[w.length()-1].find(w);
  // if we're counting ds, the string might be stored as its complement
  if((iter==wp_[w.length()-1].end())&&ds_) iter= wp_[w.length()-1].find(ReverseComplement()(w));
  if(iter==wp_[w.length()-1].end()) return 0.0;
  double sum=0.0;
  for(iVec::const_iterator i=(*iter).second.begin();i!=(*iter).second.end();i++) sum += (*i);
  return sum/static_cast<double>((*iter).second.size());
}
/////////////////////////////////////////////////////////////////////////////////////////////////
// cjPWCount Method: PosVar
//    returns the variance of position of the word, assuming normal distribution, which may not
//    be a good assumption!
double cjPWCount::PosVar(const string& w) const {
   strIVecMap::const_iterator iter= wp_[w.length()-1].find(w);
  // if we're counting ds, the string might be stored as its complement
  if((iter==wp_[w.length()-1].end())&&ds_) iter= wp_[w.length()-1].find(ReverseComplement()(w));
  if(iter==wp_[w.length()-1].end()) return 0.0;
  double sum=0.0;
  for(iVec::const_iterator i=(*iter).second.begin();i!=(*iter).second.end();i++)
    sum += (*i)*(*i);
  double mean= PosMean(w);
  return sum/static_cast<double>((*iter).second.size())-mean*mean;
} 
//////////////////////////////////////////////////////////////////////////////////////////////////
//  cjPWCount Method: ChiSqr
//     Calculates a Chi-squared value, which measures the uniformity of the distribution of the
//     word in question;  takes a number of bins or windows as an argument and returns the ChiSquare
//     value, as well as the binned or windowed counts in the reference variable bins
//     **NOTE**: bin is determined by the starting position of any word :**NOTE**
double cjPWCount::ChiSqr(const string& w,unsigned nw,dVec& bins) const {
  bins.erase(bins.begin(),bins.end()); 
  double wc= WordCount(w),dnw= nw;
  double nkMean= wc/dnw;  // mean number of counts per bin (count/number of bins)
//	WriteFloatToLog()(string("in ChiSqr,nkMean= "),nkMean);
//	WriteFloatToLog()(string("in ChiSqr,wc= "),wc);
//	WriteFloatToLog()(string("in ChiSqr,dnw= "),dnw);
  if(nkMean < 5.0) return -1.0;
  for(unsigned i=0;i<nw;i++) bins.push_back(0);
  strIVecMap::const_iterator iter= wp_[w.length()-1].find(w);
  // if the count is ds, we have to check for which word is stored
  if(iter==wp_[w.length()-1].end()&&ds_) iter= wp_[w.length()-1].find(ReverseComplement()(w));
  if(iter==wp_[w.length()-1].end()) return 0.0;
  // wSize is the size of the window for which we are binning.  noCount will count words left at
  // the end for a final incomplete window
  unsigned wSize= (maxSeqL_-w.length()+1)/nw,noCount=0;
  // j iterates through the positions
  for(iVec::const_iterator j=(*iter).second.begin();j!=(*iter).second.end();j++){
    unsigned idx=(*j)/wSize; // bin index is just position/window size (integer division, no remainder)
    if(idx<nw) bins[idx]++;
    else       noCount++;
  }
//	WriteIntToLog()(string("in ChiSqr, noCount= "),noCount);
//	WriteIntToLog()(string("in ChiSqr, minSeqL_= "),minSeqL_);
  // the chi-squared variable is  sum [(bin-mean)^2/mean]
  double sum=0.0;
  for(unsigned k=0;k<bins.size();k++){
    sum += (bins[k]-nkMean)*(bins[k]-nkMean)/nkMean;
  }
  return sum;
} 
/////////////////////////////////////////////////////////////////////////////////////////////////////
// cjPWCount Method: zChiSqr
//    same as ChiSqr, except the bins are anchored on a zero point, and bins are defined here by 
//    their length (wL) instead of their number as in the ChiSqr function
double cjPWCount::zChiSqr(const string& w,unsigned wL,int zPos,dVec& bins,double& avgBinCount) const {
  bins.erase(bins.begin(),bins.end()); 
  unsigned nw= zPos/wL + (maxSeqL_-zPos-w.length()+1)/wL; // nw = number of windows
  double wc= WordCount(w),dnw= nw; // dnw = double version of nw, wc= word count of w
  int leftOffset= (zPos-wL*(zPos/wL)); // unused positions (before first full window)
//
  for(unsigned i=0;i<nw;i++) bins.push_back(0);
  if(!WordCount(w)) return 0.0;
  // iter points into the map, to this word's vector of positions
  strIVecMap::const_iterator iter= wp_[w.length()-1].find(w);
  if(iter==wp_[w.length()-1].end()) iter= wp_[w.length()-1].find(ReverseComplement()(w));
  if(iter==wp_[w.length()-1].end()) return double(0.0);
  unsigned noCount=0;
  // loop through each position, putting it into the right bin
  for(iVec::const_iterator j=(*iter).second.begin();j!=(*iter).second.end();j++){
    if((*j)>=leftOffset) {
      unsigned idx=((*j)-leftOffset)/wL;
      if(idx<nw) //idx=nw-1;
	bins[idx]++;
      else
	noCount++;
    } else
      noCount++;
  }
  //	
  //	WriteIntToLog()(string("in zChiSqr: word: ")+w+string(" noCount="),noCount);
  avgBinCount= (wc-static_cast<double>(noCount))/dnw;  
  double sum=0.0;
  for(unsigned k=0;k<bins.size();k++){
    sum += (bins[k]-avgBinCount)*(bins[k]-avgBinCount)/avgBinCount;
  }
  return sum;
} 
////////////////////////////////////////////////////////////////////////////////////////////////
//
unsigned cjPWCount::SimProfWords(const string& w,strVec& words,dVec& aBins,unsigned nw,
				      double minR,double minChi,const strVec& fWords) const {
  words.erase(words.begin(),words.end()); aBins.erase(aBins.begin(),aBins.end());
  if(!WordCount(w)) return 0;  // check for word in the count
  //	dVec wBins;
  if(ChiSqr(w,nw,aBins)<minChi) return 0; // don't bother if the distribution isn't interesting
  words.push_back(w); bool cont=true;
  while (cont) {
    unsigned startSize= words.size();
    for(strIntMap::const_iterator it=wc_[w.length()-1].begin();it!=wc_[w.length()-1].end();it++) 
      if( (find(words.begin(),words.end(),(*it).first)==words.end())&&
	  (find(fWords.begin(),fWords.end(),(*it).first)==fWords.end())){
	dVec bins;
	if(ChiSqr((*it).first,nw,bins)>=minChi) {
	  if(Pearson()(aBins,bins)>=minR) {
	    words.push_back((*it).first);
	    for(unsigned i=0;i<aBins.size();i++) aBins[i] += bins[i];
	  }
	}
      }
    cont= (words.size()!=startSize);
  }
  return words.size();
}
///////////////////////////////////////////////////////////////////
//
void cjPWCount::GenerateRandomChi(unsigned nw,unsigned order,unsigned nTest,unsigned nSeqs,
				       unsigned seqL,dVec& chiBound) const {
  chiBound.erase(chiBound.begin(),chiBound.end());
  for(unsigned i1=0;i1<mN_;i1++) chiBound.push_back(0.0);
  for(unsigned i2=0;i2<nTest;i2++) {
    cjPWCount tempCount(mN_);
    for(unsigned j=0;j<nSeqs;j++) {
      string seq= GenerateSequence(seqL,order);
//			WriteToLog()(seq);
      tempCount.AddSequence(seq);
    }
    strVec words=tempCount.ListWords(); dVec bins;
    for(unsigned j2=0;j2<words.size();j2++) {
      double chiSq=tempCount.ChiSqr(words[j2],nw,bins);
//			WriteToLog()(string("word: ")+words[j2]);
//			WriteIntToLog()(string("count= "),tempCount.WordCount(words[j2]));
//			WriteFloatToLog()(string("mean p= "),tempCount.PosMean(words[j2]));
//			WriteFloatToLog()(string(" chi= "),chiSq);
//			for(unsigned j3=0;j3<bins.size();j3++)
//				WriteIntToLog()(string("bin: "),bins[j3]);
      if(chiSq>chiBound[words[j2].length()-1]) 
	chiBound[words[j2].length()-1]=chiSq;
    }
  }
  for(unsigned i3=0;i3<mN_;i3++) chiBound[i3] *= 1.25;
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////
// cjPWCount Method: RangeWordCount
//    returns the count of words in the range r defined by l <= r < h
unsigned cjPWCount::RangeWordCount(string& w,unsigned l,unsigned h) const {
  if(!WordCount(w)) return 0;
  // iter points into the map, to this word's vector of positions
  strIVecMap::const_iterator iter= wp_[w.length()-1].find(w);
  if((iter==wp_[w.length()-1].end())&&ds_) iter= wp_[w.length()-1].find(ReverseComplement()(w));
  if(iter==wp_[w.length()-1].end()) return 0;
  unsigned count=0;
  // loop through all positions in the vector, adding to the count if it's in the range
  for(iVec::const_iterator j=(*iter).second.begin();j!=(*iter).second.end();j++)
    if( ((*j)>=l)&&((*j)<h)) count++;
  return count;
}
///////////////////////////////////////////////////////////////////
//
void cjPWCount::KMeanCluster(unsigned k,unsigned n,unsigned wL,int zPos,double minAvgBinCount,
				  unsigned maxIt,dMat& bins,strVec& words,iVec& cID,dVec& Ssq,
				  dMat& clusterCenters,dMat& clusterPr) const {
  if(k>10) { // maximum of 10 clusters
    WriteIntToLog()(string("k for k-mean clustering reset to max (10):"),k);
    k=10;
  }
  if((n>mN_)||(n<1)) { // word size must be within the range counted 
    WriteIntToLog()(string("n for k-mean clustering reset to max:"),n);
    n=mN_;
  }
  // write info to the log file
  WriteToLog()(string("K-means parameters"));
  WriteIntToLog()(string("K-means k:"),k);	WriteIntToLog()(string("K-means n:"),n);
  WriteIntToLog()(string("window length:"),wL);	
  WriteIntToLog()(string("zero position:"),zPos);
  WriteFloatToLog()(string("MinAvgBinCount:"),minAvgBinCount);
  WriteIntToLog()(string("Max Iterations:"),maxIt);
  // clear the input arrays
  bins.erase(bins.begin(),bins.end()); words.erase(words.begin(),words.end()); 
  cID.erase(cID.begin(),cID.end()); Ssq.erase(Ssq.begin(),Ssq.end());
  // use the time as seed for the random variable
  time_t t0; time(&t0); seed= -1*t0;
  // loop through all words and get useable 
  unsigned allWords= pow(4,n); 
  for(unsigned i=0;i<allWords;i++) {
    dVec tBins; double avgCount; string word= IdxToSeq()(i,n);
    //const string& w,unsigned wL,int zPos,dVec& bins,double& avgBinCount) 
    double zChi= zChiSqr(word,wL,zPos,tBins,avgCount);
    if(avgCount >= minAvgBinCount) {
      bins.push_back(tBins); words.push_back(word); Ssq.push_back(zChi);
      // randomly assign to one of the k-clusters
      int idx= (static_cast<int>(static_cast<double>(k) * rand3()(&seed))) % k;
      cID.push_back(idx); 
    }
  }
  // now make the matrix to hold the cluster centers
  //	dMat clusterCenters;
  for(unsigned mm=0;mm<clusterCenters.size();mm++)
    clusterCenters[mm].erase(clusterCenters[mm].begin(),clusterCenters[mm].end());
  clusterCenters.erase(clusterCenters.begin(),clusterCenters.end());
  for(unsigned j=0;j<k;j++) {
    dVec temp; for(unsigned k=0;k<bins[0].size();k++) temp.push_back(0.0); 
    clusterCenters.push_back(temp);
  }
  // now do the iterative clustering
  bool cont=true; unsigned itN=0;
  while(cont) {
    // rebuild the centers;
    for(unsigned i1=0;i1<clusterCenters.size();i1++)
      for(unsigned i2=0;i2<clusterCenters[i1].size();i2++)
	clusterCenters[i1][i2]=0.0;
    for(unsigned j1=0;j1<bins.size();j1++)
      for(unsigned j2=0;j2<bins[j1].size();j2++)
	clusterCenters[cID[j1]][j2] += bins[j1][j2];
    for(unsigned m1=0;m1<clusterCenters.size();m1++) {
      double sum=0.0; 
      for(unsigned m2=0;m2<clusterCenters[m1].size();m2++) sum+=clusterCenters[m1][m2];
      for(unsigned m2=0;m2<clusterCenters[m1].size();m2++) clusterCenters[m1][m2]/=sum;
    }
    // now loop through each bin, test its distance from each cluster center, 
    // and reassign it if needed
    bool change=false;
    for(unsigned k1=0;k1<bins.size();k1++) {
      //			double minDist= EuclidDist(clusterCenters[0],bins[k1]); int min=0;
      double minDist= 1.0-Pearson()(clusterCenters[0],bins[k1]); int min=0;
      for(unsigned k2=1;k2<clusterCenters.size();k2++) {
	//				double dist=EuclidDist(clusterCenters[k2],bins[k1]);
	double dist=1.0-Pearson()(clusterCenters[k2],bins[k1]);
	if(dist<minDist) { min=k2; minDist=dist; }
      }
      if(min!=cID[k1]) {
	change=true; cID[k1]=min;
      }
    }
    if(!change || (++itN>=maxIt)) cont=false;
  }
  // now, compute the probability of assignment to each cluster, which is just 1-r(word,cluster)
  // clear the matrix first
  for(unsigned mm=0;mm<clusterPr.size();mm++)
    clusterPr[mm].erase(clusterPr[mm].begin(),clusterPr[mm].end());
  clusterPr.erase(clusterPr.begin(),clusterPr.end());
  // now loop through each word, checking against each cluster
  for(unsigned mm=0;mm<bins.size();mm++) {
    dVec tempPr;
    for(unsigned mmm=0;mmm<clusterCenters.size();mmm++)
      tempPr.push_back((1.0-Pearson()(bins[mm],clusterCenters[mmm]))/2.0);
    clusterPr.push_back(tempPr);
  }
  WriteIntToLog()(string("k-means iterations: "),itN);
  //
  return;
}
/*
////////////////////////////////////////////////////////////////////////////////////////////////
ClusterProfile::ClusterProfile(int n,const string& s,const cjWCount& wc)
:N_(n),s_(s),bc_(),wc_(),s1_(),s2_(),sz_() {
  if(N_>6) N_=6;
  int bcSize= pow(4,N_);
  bc_.erase(bc_.begin(),bc_.end()); wc_.erase(wc_.begin(),wc_.end());
  s1_.erase(s1_.begin(),s1_.end()); s2_.erase(s2_.begin(),s2_.end());
  for(unsigned i=0;i<bcSize;i++) { iVec t; wc_.push_back(t); }
  // fill the position matrix
  for(unsigned i=0;i<(s_.length()-N_+1);i++) {
    int idx= SeqToIdx()(s_.substr(i,N_));
    if(idx>=0) wc_[idx].push_back(i);
    //    if(idx>=0) {
    //      if(!wc_[idx].size()) wc_[idx].push_back(i);
    //      else if((i-wc_[idx][wc_[idx].size()-1])>=N_) wc_[idx].push_back(i);
    //    }
  }
  // now, for each word, calculate all p-values 
  for(unsigned i=0;i<bcSize;i++) {
    string w= IdxToSeq()(i,N_);
    double pWord= wc.ExpectFreq(w,1); 
    // prefix probabilities are needed for Compound Poisson calculation
    dVec prefixProbs; prefixProbs.push_back(0.0);
    for(unsigned ww=1;ww<N_;ww++) {
      string pw= w.substr(0,ww)+w; prefixProbs.push_back(wc.ExpectFreq(pw,1)/pWord);
    }
    if(wc_[i].size()<=1) {
      bc_.push_back(1.0); s1_.push_back(0); s2_.push_back(0); sz_.push_back(0); 
    } else {
      double minP=1.0; int start,stop,csize;
      for(unsigned k=wc_[i].size()-1;k>1;k--)
	for(unsigned j=0;j<(wc_[i].size()-k);j++) {
	  double x= wc_[i][j+k]-wc_[i][j]+N_;//-1;//-(k-1)*(N_-1);
	  double p= CompPoiss(w,k+1, x,prefixProbs,pWord)()/pWord;
	  if(p<0) {
	    WriteFloatToErr()(string("negative cluster p-val"),p);
	    char l[251]; std::strstream dump(l,250);
	    dump << i << ":" << IdxToSeq()(i,N_) << ", k=" << k << ", j=" << j << ", x= " << x;
	    dump << ", pWord=" << pWord;// << ", sum= " << sum;
	    dump << ", p1=" << wc_[i][j] << ", p2=" << wc_[i][j+k] << ", tc= " << wc_[i].size() << '\0';
	    WriteToErr()(string(l));
	  }
	  if(p<minP) {
	    minP= p; start= wc_[i][j]; stop= wc_[i][j+k]+N_; csize=k+1;
	  }
	}
      bc_.push_back(minP); s1_.push_back(start); s2_.push_back(stop); sz_.push_back(csize);
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////
double ClusterProfile::pearsCorr(const ClusterProfile& o) const { return Pearson()(bc_,o.bc_); }
////////////////////////////////////////////////////////////////////////////////////////////////
double ClusterProfile::spearCorr(const ClusterProfile& o) const { return Spearman()(bc_,o.bc_); }
////////////////////////////////////////////////////////////////////////////////////////////////
*/
////////////////////////////////////////////////////////////////////////////////////////////////
