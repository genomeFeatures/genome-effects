#ifndef __ARG_H
#define __ARG_H
//
#include <string>
#include <vector>
#include <map>
#include <functional>
using std::string; using std::map; using std::vector; using std::less;
//
typedef vector<string> strVec;
typedef vector<int> iVec;
typedef vector<double> dVec;
typedef map<string,string,less<string> > strStrMap;
const double KeyNotPresent= -1.0e30;
//
class cjArgs {
 protected:
  strStrMap args_;
  unsigned reqA_;
  strVec reqArgs_;
  string eName_;
 public:
  cjArgs():args_(),reqA_(),reqArgs_(),eName_() {}
  cjArgs(const cjArgs& a):args_(a.args_),reqA_(a.reqA_),reqArgs_(a.reqArgs_),eName_(a.eName_) {}
  cjArgs(int argc, char *argv[], unsigned reqA=0);
  unsigned numOptArgs() const { return args_.size(); }
  unsigned numAllArgs() const { return (args_.size()+reqArgs_.size()); }
  string Executable() const { return string(eName_); }
  strVec RequiredArgs() const { return strVec(reqArgs_); }
  unsigned nRequiredArgs() const { return reqA_; }
  bool HasRequiredArgs() const { return (reqA_==reqArgs_.size()); }
  string ArgValue(const string& k) const;
  bool KeyPresent(const string& k) const;
  bool GetFloatValue(const string& k,const string& l,double& v,double min,double max,
		     bool inclusive_low=false,bool inclusive_high = false) const;
  bool GetIntValue(const string& k,const string& l,int& v,int min,int max,
		   bool inclusive_low=false,bool inclusive_high = false) const;
  bool GetUnsignedValue(const string& k,const string& l,unsigned& v,unsigned min,unsigned max,
			bool inclusive_low=false,bool inclusive_high = false) const;
  bool GetString(const string& k,const string& l,string& v) const;
  bool GetFilestring(const string& k,const string& l,string& v) const;
  bool GetIntVec(const string& k,const string& l,iVec& iv) const;
  bool GetFloatVec(const string& k,const string& l,dVec& iv) const;
};
//
#endif
