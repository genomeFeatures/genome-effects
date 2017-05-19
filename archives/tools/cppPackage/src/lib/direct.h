////////////////////////////////////////////////////////////////////////////////////////////
#ifndef __direct_h
#define __direct_h
#include "util.h"
#include <dirent.h>
#include <string>
#include <vector>
#include <functional>
#include <map>
////////////////////////////////////////////////////////////////////////////////////////////
using std::string; using std::vector; using std::map;
typedef vector<string> strVec;
typedef map<string,dirent *,std::less<string> > strFptrMap;
////////////////////////////////////////////////////////////////////////////////////////////
class directory {
 private:
  string path_;
  strFptrMap fMap_;
 public:
  directory();
  directory(const string& path);
  directory(const directory& d);
  //
  bool hasEntry(const string& pat, bool noHidden=true) const;
  bool hasFile(const string& pat, bool noHidden=true) const;
  bool hasSubDirectory(const string& pat, bool noHidden=true) const;
  //
  strVec SubDirectoryList(const string& pat=string(),bool noHidden=true) const;
  strVec FileList(const string& pat=string(),bool noHidden=true) const;
  strVec FullList(bool noHidden=true) const;
  //
  string Path() const { return path_; }
};
#endif
