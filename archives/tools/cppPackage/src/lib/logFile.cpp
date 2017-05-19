// written and debugged by Joel Graber, Associate Staff Scientist at the Jackson Laboratory
// All rights reserved, 1999-2003.
//
// logFile.cpp - unary function that writes a string to the log file, with a timestamp
// creating the log file if necessary.
#include "../lib/logFile.h"
#include <time.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
//////////////////////////////////////////////////////////////////////////////////////////////
std::string logFileName("default.log");
std::string errFileName("default.err");
std::string LogFileName() { return logFileName; }
std::string ErrFileName() { return errFileName; }
//////////////////////////////////////////////////////////////////////////////////////////////
// Global variable WroteError_- indicates error file activity other than initiation
// Global variable LogsOff_   - a way of shutting off logging altogether for performance
bool WroteError_=false,LogsOff_=false;
bool WroteError() { return WroteError_; } //
void LogsOff() { LogsOff_=true; return; }
//////////////////////////////////////////////////////////////////////////////////////////////
//
bool LogSetup::operator()(int argc,char *argv[]) const {
  if(LogsOff_) return true;
  strVec args; args.reserve(argc); for(int i=0;i<argc;++i) args.push_back(std::string(argv[i]));
  bool RunOK= WriteRunLog()(args);
  if(!RunOK) {
    chdir("/tmp/logs/");
    RunOK= WriteRunLog()(args);
    if(!RunOK) return false;
  }
  char line[501]; getcwd(line,500); std::string oldPath(line);
  if(oldPath[oldPath.length()-1]!='/') oldPath += '/';
  std::string arg0(argv[0]);  int s=arg0.rfind('/'); 
  std::string prName(arg0.substr(s+1,arg0.length()-s-1));
  bool LogOK= InitLog()(oldPath+prName)==1;
  bool ErrOK= InitErrLog()(oldPath+prName)==1;
  if(!LogOK || !ErrOK) {
    oldPath= std::string("/tmp/logs/");
    LogOK= InitLog()(oldPath+prName)==1;
    ErrOK= InitErrLog()(oldPath+prName)==1;
  }
  return LogOK&&ErrOK;
}  
//////////////////////////////////////////////////////////////////////////////////////////////
//
time_t t0,t1; // global variables to hold logfile initiation time.
int InitLog::operator()(const std::string& fn,bool overwrite) const {
  if(LogsOff_) return 1;
  logFileName=fn+std::string(".log");  
  std::ofstream logFile; 
  if(overwrite) logFile.open(logFileName.c_str(),std::ios::out);
  else          logFile.open(logFileName.c_str(),std::ios::app);
  if(logFile.fail()) return(2);
  time(&t0); t1=t0;
  if(overwrite) {
    logFile << ctime(&t0) << "Program start and initialization of log" << std::endl;
    logFile << "code written by Joel Graber (jhgraber@jax.org; 207-288-6782)" << std::endl;
    logFile << "Associate Staff Scientist at The Jackson Laboratory" << std::endl;
  }
  logFile.close();
  return(1);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//  timestamp write to log
int WriteToLogTS::operator()(const std::string& s) const {
  if(LogsOff_) return 1;
  std::ofstream logFile; logFile.open(logFileName.c_str(),std::ios::app);
  if(logFile.fail()) logFile.open(logFileName.c_str(),std::ios::out);
  if(logFile.fail()) return(2);
  time_t t; time(&t);
  logFile << s.c_str() << ":" << ctime(&t);
  logFile << (t-t1) << " s elapsed since last TS, ";
  logFile << (t-t0) << " s elapsed since startup." << std::endl;
  logFile.close(); t1=t;
  return(1);
}//////////////////////////////////////////////////////////////////////////////////////////////
//  timestamp write to log
int CloseLogTS::operator()(const std::string& s) const {
  if(LogsOff_) return 1;
  std::ofstream logFile; logFile.open(logFileName.c_str(),std::ios::app);
  if(logFile.fail()) logFile.open(logFileName.c_str(),std::ios::out);
  if(logFile.fail()) return(2);
  time_t t; time(&t);
  logFile << s.c_str() << ":" << ctime(&t);
  logFile << (t-t1) << " s elapsed since last TS, ";
  logFile << (t-t0) << " s elapsed since startup." << std::endl;
  logFile.close(); t1=t;
  return(1);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
int WriteToLog::operator()(const std::string& s) const {
  if(LogsOff_) return 1;
  std::ofstream logFile; logFile.open(logFileName.c_str(),std::ios::app);
  if(logFile.fail()) logFile.open(logFileName.c_str(),std::ios::out);
  if(logFile.fail()) return(2);
  logFile << s.c_str() << std::endl;
  logFile.close();
  return(1);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
int WriteIntToLog::operator()(const std::string& s,int i) const {
  if(LogsOff_) return 1;
  std::ofstream logFile; logFile.open(logFileName.c_str(),std::ios::app);
  if(logFile.fail()) logFile.open(logFileName.c_str(),std::ios::out);
  if(logFile.fail()) return(2);
  logFile << s.c_str() << " " << i << std::endl;
  logFile.close();
  return(1);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
int WriteFloatToLog::operator()(const std::string& s,float f) const {
  if(LogsOff_) return 1;
  std::ofstream logFile; logFile.open(logFileName.c_str(),std::ios::app);
  if(logFile.fail()) logFile.open(logFileName.c_str(),std::ios::out);
  if(logFile.fail()) return(2);
  logFile << s.c_str() << " " << f << std::endl;
  logFile.close();
  return(1);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
time_t e0,e1; // global variables to hold logfile initiation time.
int InitErrLog::operator()(const std::string& fn,bool overwrite) const {
  if(LogsOff_) return 1;
  errFileName=fn+std::string(".err");
  std::ofstream logFile; 
  if(overwrite) logFile.open(errFileName.c_str(),std::ios::out);
  else          logFile.open(errFileName.c_str(),std::ios::app);
  if(logFile.fail()) return(2);
  time(&e0); e1=e0;
  if(overwrite) {
    logFile << ctime(&t0) << "Program start and initialization of log" << std::endl;
    logFile << "code written by Joel Graber (jhgraber@jax.org; 207-288-6782)" << std::endl;
    logFile << "Associate Staff Scientist at The Jackson Laboratory" << std::endl;
  }
  logFile.close();
  return(1);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//  timestamp write to log
int WriteToErrTS::operator()(const std::string& s) const {
  if(LogsOff_) return 1;
  std::ofstream logFile; logFile.open(errFileName.c_str(),std::ios::app);
  if(logFile.fail()) logFile.open(errFileName.c_str(),std::ios::out);
  if(logFile.fail()) return(2);
  time_t e; time(&e);
  logFile << s.c_str() << ":" << ctime(&e);
  logFile << (e-e1) << " s elapsed since last TS, ";
  logFile << (e-e0) << " s elapsed since startup." << std::endl;
  logFile.close(); e1=e;
  WroteError_=true;
  return(1);
}//////////////////////////////////////////////////////////////////////////////////////////////
//  timestamp write to log
int CloseErrTS::operator()(const std::string& s) const {
  if(LogsOff_) return 1;
  std::ofstream logFile; logFile.open(errFileName.c_str(),std::ios::app);
  if(logFile.fail()) logFile.open(errFileName.c_str(),std::ios::out);
  if(logFile.fail()) return(2);
  time_t e; time(&e);
  logFile << s.c_str() << ":" << ctime(&e);
  logFile << (e-e1) << " s elapsed since last TS, ";
  logFile << (e-e0) << " s elapsed since startup." << std::endl;
  logFile.close(); e1=e;
  return(1);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
int WriteToErr::operator()(const std::string& s) const {
  if(LogsOff_) return 1;
  std::ofstream logFile; logFile.open(errFileName.c_str(),std::ios::app);
  if(logFile.fail()) logFile.open(errFileName.c_str(),std::ios::out);
  if(logFile.fail()) return(2);
  logFile << s.c_str() << std::endl;
  logFile.close();
  WroteError_=true;
  return(1);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
int WriteIntToErr::operator()(const std::string& s,int i) const {
  if(LogsOff_) return 1;
  std::ofstream logFile; logFile.open(errFileName.c_str(),std::ios::app);
  if(logFile.fail()) logFile.open(errFileName.c_str(),std::ios::out);
  if(logFile.fail()) return(2);
  logFile << s.c_str() << " " << i << std::endl;
  logFile.close();
  WroteError_=true;
  return(1);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
int WriteFloatToErr::operator()(const std::string& s,float f) const {
  if(LogsOff_) return 1;
  std::ofstream logFile; logFile.open(errFileName.c_str(),std::ios::app);
  if(logFile.fail()) logFile.open(errFileName.c_str(),std::ios::out);
  if(logFile.fail()) return(2);
  logFile << s.c_str() << " " << f << std::endl;
  logFile.close();
  WroteError_=true;
  return(1);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
int WriteToLogs::operator()(const std::string& s) const {
  if(LogsOff_) return 1;
  std::ofstream logFile; logFile.open(logFileName.c_str(),std::ios::app);
  if(logFile.fail()) logFile.open(logFileName.c_str(),std::ios::out);
  if(logFile.fail()) return(2);
  logFile << s.c_str() << std::endl;
  logFile.close();
//
  logFile.open(errFileName.c_str(),std::ios::app);
  if(logFile.fail()) logFile.open(errFileName.c_str(),std::ios::out);
  if(logFile.fail()) return(2);
  logFile << s.c_str() << std::endl;
  logFile.close();
  WroteError_=true;
  //
  return(1);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
int WriteIntToLogs::operator()(const std::string& s,int i) const {
  if(LogsOff_) return 1;
  std::ofstream logFile; logFile.open(logFileName.c_str(),std::ios::app);
  if(logFile.fail()) logFile.open(logFileName.c_str(),std::ios::out);
  if(logFile.fail()) return(2);
  logFile << s.c_str() << " " << i << std::endl;
  logFile.close();
//
  logFile.open(errFileName.c_str(),std::ios::app);
  if(logFile.fail()) logFile.open(errFileName.c_str(),std::ios::out);
  if(logFile.fail()) return(2);
  logFile << s.c_str() << " " << i << std::endl;
  logFile.close();
  WroteError_=true;
  return(1);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
int WriteFloatToLogs::operator()(const std::string& s,float f) const {
  if(LogsOff_) return 1;
  std::ofstream logFile; logFile.open(logFileName.c_str(),std::ios::app);
  if(logFile.fail()) logFile.open(logFileName.c_str(),std::ios::out);
  if(logFile.fail()) return(2);
  logFile << s.c_str() << " " << f << std::endl;
  logFile.close();
  //
  logFile.open(errFileName.c_str(),std::ios::app);
  if(logFile.fail()) logFile.open(errFileName.c_str(),std::ios::out);
  if(logFile.fail()) return(2);
  logFile << s.c_str() << " " << f << std::endl;
  logFile.close();
  WroteError_=true;
  return(1);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
bool WriteRunLog::operator()(const strVec& a) const {
  if(LogsOff_) return true;
  time_t t; time(&t);
  std::ofstream runLog;  
  runLog.open("run.log",std::ios::app);
  if(runLog.fail()) {
    runLog.open("run.log",std::ios::out);
  }
  if(runLog.fail()) {
    return false;
  }
  for(unsigned i=0;i<a.size();i++)
    if(!i) runLog << a[i];
    else  runLog << ' ' << a[i];
  runLog << ": " << ctime(&t);
  runLog.close();
  return true;
}
  
