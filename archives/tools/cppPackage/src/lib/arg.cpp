#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <strstream>
#include <iostream>
#include "../lib/arg.h"
#include "../lib/logFile.h"
#include "../lib/util.h"
/////////////////////////////////////////////////////////////////////////////////////
//
cjArgs::cjArgs(int argc, char *argv[], unsigned reqA):reqA_(reqA),args_(),reqArgs_(),eName_() {
  int dbcount=0; 
  //  std::cout << "argc= " << argc << ", reqA= " << reqA << std::endl;
  WriteIntToLog()(std::string("cjArgs"),++dbcount);
  eName_= std::string(argv[0]);
  if(argc<(reqA+1)) return;  // program name + required arguments = minimum acceptable args
  for(unsigned i=0;i<reqA;i++) 
    reqArgs_.push_back(std::string(argv[i+1])); // put the required arguments into the vector
  // create two arrays- a string array of all optional arguments
  //  and an int array of flags, set to positive if the corresponding optArg starts with '-'
  //  unless the next character is a digit (indicating a negative number)
  strVec optArgs; optArgs.reserve(argc-1-reqA);
  iVec flag(argc-1-reqA,0);
  for(unsigned j=reqA+1;j<argc;++j) {
    optArgs.push_back(std::string(argv[j]));
    if(optArgs.back()[0]=='-') {
      if(optArgs.back().size()>1 && !isdigit(optArgs.back()[1])) flag[j-reqA-1]=1;
    }
  }
  // now process the optional arguments, putting pairs into the string-string map args_
  unsigned j=0;
  while(j<optArgs.size()) {
    if(!flag[j]) {
      WriteToLogs()(std::string("invalid option: ")+std::string(argv[j]));
      ++j;
    } else if(optArgs[j].size()==1) {
      WriteToLogs()(std::string("bare dash in arguments"));
      ++j;
    } else { // optional args can be one of two forms -k v or -k=v (no spaces allowed
      //check for '=' character, 
      strVec f=split()(optArgs[j],'=');
      std::string key,val;
      if(f.size()==1) { // no '=' found
	key= optArgs[j].substr(1,optArgs[j].size()-1);
	++j;
	if(j<optArgs.size() && !flag[j])
	  val= optArgs[j++];
      } else {
	key= f[0];
	val= f[1];
	++j;
      }
      /*  unsigned j=reqA+1; // check for and process optional arguments
  while(j<argc) {
    if(argv[j][0] != '-') { // all optional arguments begin with a '-' character
      std::cout << "dash " << (++dbcount) << std::endl;
      WriteToErr()(std::string("invalid option: ")+std::string(argv[j]));
      WriteToLog()(std::string("invalid option: ")+std::string(argv[j]));
      ++j;
    } else {
      std::cout << "other " << (++dbcount) << std::endl;
      std::string full(argv[j]),k,v;  // two allowed ways of key-val args
      strVec f= split(full,'=');
      if(f.size()==1) { // no '=' found, get value from the next field 
	k= full.substr(1,full.length()-1);
	 
      unsigned split=full.find("=");  // -k=v or -k v
      if(split<full.size()) { // = found key and val are together
	k= full.substr(1,split-1);
	v= full.substr(split+1,full.length()-split-1);
	std::cout << "split " << (++dbcount) << std::endl;
	++j;
      } else {
	k= full.substr(1,full.length()-1);
	if(j<argc) {
	  if(argv[j][0]!='-')
	    v= std::string(argv[j++]);
	}
	std::cout << "nosplit " << (++dbcount) << std::endl;
	} */
      if(args_.count(key)) {
	//std::cout << "count " << (++dbcount) << std::endl;
	WriteToErr()(std::string("duplicate optional arg ignored: ")+key);
	WriteToLog()(std::string("duplicate optional arg ignored: ")+key);
      } else {
	//std::cout << "nocount " << (++dbcount) << ", k=" << key << ", j= " << j << std::endl;
	WriteToLog()(std::string("command line k-v pair: ")+key+std::string(":")+val);
	args_[key]= val;
      }
    } //else
  } // end of while loop )
  return;
} // cjArgs::cjArgs(int argc, char *argv[], unsigned reqA)
/////////////////////////////////////////////////////////////////////////////////////
//
bool cjArgs::KeyPresent(const std::string& key) const {
  return (args_.find(key)!=args_.end());
}// bool cjArgs::KeyPresent(const std::string& key) 
/////////////////////////////////////////////////////////////////////////////////////
//
std::string cjArgs::ArgValue(const std::string& key) const {
  if(KeyPresent(key)) {
    WriteToLog()(std::string("Argument: ")+key+std::string(", value: ")+(*args_.find(key)).second);
    return (*args_.find(key)).second;
  }
  return std::string();
} // std::string cjArgs::ArgValue(const std::string& key)
/////////////////////////////////////////////////////////////////////////////////////
bool cjArgs::GetFloatValue(const std::string& k,const std::string& l,double& v,double n,double x,
			   bool inclusive_low,bool inclusive_high) const {
  if(!KeyPresent(k)) {
    WriteFloatToLog()(std::string("default value of ")+l+std::string(" = "),v);
    return false;
  }
  double t= atof((*args_.find(k)).second.c_str());
  if((t<n)||(t>x)||(!inclusive_low&&(t==n))||(!inclusive_high&&(t==x))) {
    char line[301]; std::strstream ll(line,300);
    ll << l.c_str() << " invalid float input: " << t << ", reset to " << v;
    if(inclusive_low&&inclusive_high) ll << " [range: " << n << " < v < " << x << "]";
    else if(inclusive_low)            ll << " [range: " << n << " <= v < " << x << "]";
    else if(inclusive_high)           ll << " [range: " << n << " < v <= " << x << "]";
    else                              ll << " [range: " << n << " <= v <= " << x << "]";
    ll << '\0';
    WriteToLogs()(std::string(line));
    return false;
  }
  v=t;
  WriteFloatToLog()(std::string("double: ")+l+std::string(" = "),v);
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////
bool cjArgs::GetIntValue(const std::string& k,const std::string& l,int& v,int n,int x,
			 bool inclusive_low,bool inclusive_high) const {
  if(!KeyPresent(k)) {
    WriteIntToLog()(std::string("default value of ")+l+std::string(" = "),v);
    return false;
  }
  int t= atoi((*args_.find(k)).second.c_str());
  if((t<n)||(t>x)||(!inclusive_low&&(t==n))||(!inclusive_high&&(t==x))) {
    char line[301]; std::strstream ll(line,300);
    ll << l.c_str() << " invalid int input: " << t << ", reset to " << v;
    if(inclusive_low&&inclusive_high) ll << " [range: " << n << " < v < " << x << "]";
    else if(inclusive_low)            ll << " [range: " << n << " <= v < " << x << "]";
    else if(inclusive_high)           ll << " [range: " << n << " < v <= " << x << "]";
    else                              ll << " [range: " << n << " <= v <= " << x << "]";
    ll << '\0';
    WriteToLogs()(std::string(line));
    return false;
  }
  v=t;
  WriteIntToLog()(std::string("int: ")+l+std::string(" = "),v);
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////
bool cjArgs::GetUnsignedValue(const std::string& k,const std::string& l,unsigned& v,unsigned n,unsigned x,
			      bool inclusive_low,bool inclusive_high) const {
  if(!KeyPresent(k)) {
    WriteIntToLog()(std::string("default value of ")+l+std::string(" = "),v);
    return false;
  }
  unsigned t= atoi((*args_.find(k)).second.c_str());
  if((t<n)||(t>x)||(inclusive_low&&(t==n))||(inclusive_high&&(t==x))) {
    char line[301]; std::strstream ll(line,300);
    ll << l.c_str() << " invalid unsigned input: " << t << ", reset to " << v;
    if(inclusive_low&&inclusive_high) ll << " [range: " << n << " < v < " << x << "]";
    else if(inclusive_low)            ll << " [range: " << n << " <= v < " << x << "]";
    else if(inclusive_high)           ll << " [range: " << n << " < v <= " << x << "]";
    else                              ll << " [range: " << n << " <= v <= " << x << "]";
    ll << '\0';
    WriteToLogs()(std::string(line));
    return false;
  }
  v=t;
  WriteIntToLog()(std::string("unsigned: ")+l+std::string(" = "),v);
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////
bool cjArgs::GetString(const std::string& k,const std::string& l,std::string& v) const {
  if(!KeyPresent(k)) {
    WriteToLog()(std::string("default value of ")+l+std::string(" = ")+v);
    return false;
  }
  v=(*args_.find(k)).second;
  WriteToLog()(std::string("std::string: ")+l+(" = ")+v);
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////
bool cjArgs::GetFilestring(const std::string& k,const std::string& l,std::string& v) const {
  if(!KeyPresent(k)) {
    WriteToLog()(std::string("default value of ")+l+std::string(" = ")+v);
    return false;
  }
  std::string fnCheck= MakeFilenameSafe()((*args_.find(k)).second);
  if(fnCheck.length()) {
    v=fnCheck;
    WriteToLog()(std::string("std::string: ")+l+(" = ")+v);
  } else 
    WriteToLogs()(std::string("proposed filestd::string returned empty std::string, retaining default: ")+v);  
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////
bool cjArgs::GetIntVec(const std::string& k,const std::string& l,iVec& iv) const {
  if(!KeyPresent(k)) {
    std::string outvec;
    for(int i=0;i<iv.size();++i) outvec += std::string(" ")+IntToString()(iv[i]);
    WriteToLog()(std::string("default value of ")+l+std::string(" = ")+outvec);
    return false;
  }
  strVec f=split()((*args_.find(k)).second,':');
  if(f.empty()) {
    WriteToLog()(std::string("empty list entered for argument key ")+k);
    return false;
  }
  iv.clear(); iv.reserve(f.size());
  for(int i=0;i<f.size();++i) iv.push_back(atoi(f[i].c_str()));
  std::string outvec;
  for(int i=0;i<iv.size();++i) outvec += std::string(" ")+IntToString()(iv[i]);
  WriteToLog()(std::string("new value of ")+l+std::string(" = ")+outvec);
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////
bool cjArgs::GetFloatVec(const std::string& k,const std::string& l,dVec& iv) const {
  if(!KeyPresent(k)) {
    std::string outvec;
    for(int i=0;i<iv.size();++i) outvec += std::string(" ")+DoubleToString()(iv[i]);
    WriteToLog()(std::string("default value of ")+l+std::string(" = ")+outvec);
    return false;
  }
  strVec f=split()((*args_.find(k)).second,':');
  if(f.empty()) {
    WriteToLog()(std::string("empty list entered for argument key ")+k);
    return false;
  }
  iv.clear(); iv.reserve(f.size());
  for(int i=0;i<f.size();++i) iv.push_back(atof(f[i].c_str()));
  std::string outvec;
  for(int i=0;i<iv.size();++i) outvec += std::string(" ")+DoubleToString()(iv[i]);
  WriteToLog()(std::string("new value of ")+l+string(" = ")+outvec);
  return true;
}
