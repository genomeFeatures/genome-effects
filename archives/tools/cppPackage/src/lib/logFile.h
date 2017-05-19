// written and debugged by Joel Graber, Senior Research Associate at the
// Center for Advanced Biotechnology, Boston University
// All rights reserved by The Trustees of Boston University, 1999.
//
#ifndef __LOGFILE_H
#define __LOGFILE_H
//
#ifndef _STRING_
#include <string>
#endif
#ifndef _FUNCTIONAL_
#include <functional>
#endif
#include <vector>
typedef std::vector<std::string> strVec;
//////////////////////////////////////////////////////////////////////////////////////////////
class LogSetup:public std::binary_function<int,char**,bool> {
 public:
  bool operator()(int argc,char *argv[]) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class InitLog:public std::binary_function<std::string, bool, int> {
 public:
  int operator()(const std::string& fn,bool overwrite=true) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class InitErrLog:public std::binary_function<std::string, bool, int> {
 public:
  int operator()(const std::string& fn,bool overwrite=true) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteToLogTS:public std::unary_function<std::string,int> {
 public:
  int operator()(const std::string& s) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class CloseLogTS:public std::unary_function<std::string,int> {
 public:
  int operator()(const std::string& s) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteToLog:public std::unary_function<std::string,int> {
 public:
  int operator()(const std::string& s) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteIntToLog:public std::binary_function<std::string,int,int> {
 public:
  int operator()(const std::string& s,int i) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteFloatToLog:public std::binary_function<std::string,float,int> {
 public:
  int operator()(const std::string& s,float f) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteToErrTS:public std::unary_function<std::string,int> {
 public:
  int operator()(const std::string& s) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class CloseErrTS:public std::unary_function<std::string,int> {
 public:
  int operator()(const std::string& s) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteToErr:public std::unary_function<std::string,int> {
 public:
  int operator()(const std::string& s) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteIntToErr:public std::binary_function<std::string,int,int> {
 public:
  int operator()(const std::string& s,int i) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteFloatToErr:public std::binary_function<std::string,float,int> {
 public:
  int operator()(const std::string& s,float f) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteToLogs:public std::unary_function<std::string,int> {
 public:
  int operator()(const std::string& s) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteIntToLogs:public std::binary_function<std::string,int,int> {
 public:
  int operator()(const std::string& s,int i) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteFloatToLogs:public std::binary_function<std::string,float,int> {
 public:
  int operator()(const std::string& s,float f) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteRunLog:public std::unary_function<strVec,bool> {
 public:
  bool operator()(const strVec& a) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
bool WroteError();
void LogsOff();
std::string LogFileName();
std::string ErrFileName();
//////////////////////////////////////////////////////////////////////////////////////////////
#endif  // __LOGFILE_H
