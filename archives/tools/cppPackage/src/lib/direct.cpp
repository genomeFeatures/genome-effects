///////////////////////////////////////////////////////////////////////////////////
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include "../lib/direct.h"
#include "../lib/logFile.h"
///////////////////////////////////////////////////////////////////////////////////
bool match(const string& p,const string& s) {
  return (p==s);
}
///////////////////////////////////////////////////////////////////////////////////
directory::directory():path_(),fMap_() { 
  char line[501]; getcwd(line,500);
  path_= string(line);
  if(*(path_.rbegin())!='/') path_ += '/';
  DIR *dir= opendir(path_.c_str());
  if(!dir) {
    WriteToLogs()(string("error opening default directory: ")+path_);
    path_.clear();
    return;
  }
  struct dirent *entry;
  while((entry= readdir(dir))!=NULL) 
    fMap_[string(entry->d_name)]= entry;
  closedir(dir);
}
///////////////////////////////////////////////////////////////////////////////////
directory::directory(const string& path):path_(path),fMap_() { 
  if(path_.empty()) { char line[501]; getcwd(line,500); path_= string(line); }
  if(*(path_.rbegin())!='/') path_ += '/';
  DIR *dir= opendir(path_.c_str());
  if(!dir) {
    WriteToLogs()(string("error opening directory: ")+path_);
    path_.clear();
    return;
  }
  struct dirent *entry;
  while((entry= readdir(dir))!=NULL) 
    fMap_[string(entry->d_name)]= entry;
  closedir(dir);
}
///////////////////////////////////////////////////////////////////////////////////
directory::directory(const directory& d):path_(d.path_),fMap_(d.fMap_) {}
///////////////////////////////////////////////////////////////////////////////////
bool directory::hasEntry(const string& pat,bool noHidden) const {
  if(path_.empty()) return false;
  for(strFptrMap::const_iterator i=fMap_.begin();i!=fMap_.end();++i)
    if(match(pat,(*i).first) && (*((*i).first.begin())!='.'|| !noHidden)) 
      return true;
  return false;
}
///////////////////////////////////////////////////////////////////////////////////
bool directory::hasFile(const string& pat,bool noHidden) const {
  if(path_.empty()) return false;
  for(strFptrMap::const_iterator i=fMap_.begin();i!=fMap_.end();++i)
    if((*i).second->d_type==DT_REG && match(pat,(*i).first) && (*((*i).first.begin())!='.'|| !noHidden)) 
      return true;
  return false;
}
///////////////////////////////////////////////////////////////////////////////////
bool directory::hasSubDirectory(const string& pat,bool noHidden) const {
  if(path_.empty()) return false;
  for(strFptrMap::const_iterator i=fMap_.begin();i!=fMap_.end();++i)
    if((*i).second->d_type==DT_DIR && match(pat,(*i).first) && (*((*i).first.begin())!='.'|| !noHidden)) 
      return true;
  return false;
}
///////////////////////////////////////////////////////////////////////////////////
strVec directory::FileList(const string& pat,bool noHidden) const {
  strVec fl;
  if(!path_.empty()) 
    for(strFptrMap::const_iterator i=fMap_.begin();i!=fMap_.end();++i)
      if((*i).second->d_type==DT_REG && (pat.empty() || match(pat,(*i).first)) && (*((*i).first.begin())!='.'|| !noHidden))
	fl.push_back((*i).first);
  return fl;
}
///////////////////////////////////////////////////////////////////////////////////
strVec directory::SubDirectoryList(const string& pat,bool noHidden) const {
  strVec sd; WriteIntToLog()(string("pat= ")+pat+string(", l="),pat.length());
  if(!path_.empty()) 
    for(strFptrMap::const_iterator i=fMap_.begin();i!=fMap_.end();++i) {
      if((*i).second->d_type==DT_DIR && (pat.empty() || match(pat,(*i).first)) && (*((*i).first.begin())!='.'|| !noHidden))
	sd.push_back((*i).first);
    }
  return sd;
}
///////////////////////////////////////////////////////////////////////////////////
strVec directory::FullList(bool noHidden) const {
  strVec sd; 
  if(!path_.empty()) 
    for(strFptrMap::const_iterator i=fMap_.begin();i!=fMap_.end();++i) 
      if(!noHidden || *((*i).first.begin()) !='.')
	sd.push_back((*i).first);
  return sd;
}
///////////////////////////////////////////////////////////////////////////////////
