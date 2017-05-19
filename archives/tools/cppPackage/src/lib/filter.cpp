#include "../lib/filter.h"
#include "../lib/basicSeq.h"
#include <ctype.h>
#include <map>
#include <vector>
typedef std::map<std::string,unsigned,std::less<std::string> > strUMap;
typedef std::vector<unsigned> uVec;
//////////////////////////////////////////////////////////////////////////////////////////////
// global variables for dust routines- all dust routines downloaded from ncbi ftp site
static int word = 3; 
static int window = 64; 
static int window2 = 32; 
static int level = 20;
static int mv, iv, jv;
static char mask= 'N';
/////////////////////////////
bool SetDustLevel::operator()(int value) const { 
	if(level<0) return false;
	level = value; return true;
}
bool SetDustWindow::operator()(int value) const { 
	if(value<0) return false;
	window = value; window2 = window / 2; return true;
}
bool SetDustWord::operator()(int value) const { 
	if(level<0) return false;
	word = value; return true;
}
bool SetDustMask::operator()(char m) const {
  if(((m>='a')&&(m<='z')) || ((m>='A')&&(m<='Z'))) {
    mask=m; return true;
  }
  return false;
}
/////////////////////////////
/////////////////////////////
static void wo1(int len, char *s, int ivv) {
	int i, ii, j, v, t, n, n1, sum;
	static int counts[32*32*32];
	static int iis[32*32*32];
	int js, nis;

	n = 32 * 32 * 32;
	n1 = n - 1;
	nis = 0;
	i = 0;
	ii = 0;
	sum = 0;
	v = 0;
	for (j=0; j < len; j++, s++) {
		ii <<= 5;
		if (isalpha(*s)) {
			if (islower(*s)) {
				ii |= *s - 'a';
			} else {
				ii |= *s - 'A';
			}
		} else {
			i = 0;
			continue;
		}
		ii &= n1;
		i++;
		if (i >= word) {
			for (js=0; js < nis && iis[js] != ii; js++) ;
			if (js == nis) {
				iis[nis] = ii;
				counts[ii] = 0;
				nis++;
			}
			if ((t = counts[ii]) > 0) {
				sum += t;
				v = 10 * sum / j;
				if (mv < v) {
					mv = v;
					iv = ivv;
					jv = j;
				}
			}
			counts[ii]++;
		}
	}
}
/////////////////////////////
static int wo(int len,char *s,int *beg,int *end){
	int i, l1;

	l1 = len - word + 1;
	if (l1 < 0) {
		*beg = 0;
		*end = len - 1;
		return 0;
	}
	mv = 0;
	iv = 0;
	jv = 0;
	for (i=0; i < l1; i++) {
		wo1(len-i, s+i, i);
	}
	*beg = iv;
	*end = iv + jv;
	return mv;
}
/////////////////////////////
void dust(int len,char *s){
	int i, j, l, from, to, a, b, v;

	from = 0;
	to = -1;
	for (i=0; i < len; i += window2) {
		from -= window2;
		to -= window2;
		l = (len > i+window) ? window : len-i;
		v = wo(l, s+i, &a, &b);
		for (j = from; j <= to; j++) {
			if(DNAbases.find(s[i+j])<4 || RNAbases.find(s[i+j]<4)) s[i+j]= mask;
		}
		if (v > level) {
			for (j = a; j <= b && j < window2; j++) {
				if(DNAbases.find(s[i+j])<4 || RNAbases.find(s[i+j]<4)) s[i+j] = mask;
			}
			from = j;
			to = b;
		} else {
			from = 0;
			to = -1;
		}
	}
}
////////////////////////////////////////
std::string Dust::operator()(const std::string& s) const {
	char *ss= new char [s.length()+1];
	for(unsigned i=0;i<s.length();i++) ss[i]=s[i];
	ss[s.length()]='\0';
	dust(s.length(),ss);
	std::string ret; ret.reserve(s.length());
	ret= std::string(ss);
	delete[] ss;
	return ret;
}
//////////////////////////////////////////////////////////////////////////////////////////////
// global variables for dynamic programming filtering
static int dpWindow=64;
static int dpWord=3;
static int dpTh=1;
static int dpCut=10;
static dpFilterType dpType=sum; 
//////////////////////////////////////////////////////////////////////////////////////////////
// 
std::string DPFilter::operator()(const std::string& s) const {
	int h=0,last0=0,hMax=0,ihMax;
	uVec start,stop; strUMap wCount,hwCount; 
	unsigned seqStop=s.length()-dpWord;
	for(unsigned i=0;i<seqStop;i++) {
		std::string wR(s.substr(i,dpWord));
		int idxR=SeqToIdx()(wR);
		if(idxR>=0) {
			if(!wCount.count(wR)) wCount[wR]=0;
			wCount[wR]++; 
		}
		if(i>=dpWindow) {
			std::string wL(s.substr(i-dpWindow,dpWord));
			int idxL=SeqToIdx()(s.substr(i-dpWindow,dpWord));
			if(idxL>=0) {
				wCount[wL]--;
				if(wCount[wL]>dpTh) h -= (wCount[wL]-dpTh);
			}
			if(h<0) h=0;
		}
		if((h+wCount[wR]-1)>dpTh) {
			h += wCount[wR]-dpTh-1;
			if(h>=hMax) { hMax=h; ihMax=i; hwCount=wCount; }
		} else 			
			h=0; 
		if((hMax>dpCut)&&(!h || (i==(seqStop-1)))) {
			//find the first duplicate word in the current window
			int tempStart=0;
			if(ihMax>dpWindow) tempStart=ihMax-dpWindow;
			while((hwCount[s.substr(tempStart,dpWord)]<=dpTh)&&(tempStart<last0))tempStart++;
//			start.push_back(last0+1);
			start.push_back(tempStart);
			stop.push_back(ihMax+dpWord);
			hMax=0; 
		}		
		if(!h) { last0=i; hMax=0; }
	}
	std::string ss(s);
	for(unsigned k=0;k<start.size();k++)
		for(unsigned p=start[k];p<stop[k];p++)
			ss[p]=mask;
	return std::string(ss);
}
//////////////////////////////////////////////////////////////////////////////////////////////
bool SetDPWindow::operator()(int value) const { 
	if(value<0) return false;
	dpWindow = value; return true;
}
bool SetDPWord::operator()(int value) const { 
	if(value<0) return false;
	dpWord = value; return true;
}
bool SetDPThreshold::operator()(int value) const { 
	if(value<0) return false;
	dpTh = value; return true;
}
bool SetDPCutoff::operator()(int value) const { 
	if(value<0) return false;
	dpCut = value; return true;
}

