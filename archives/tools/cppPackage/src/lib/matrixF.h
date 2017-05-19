/////////////////// matrixF.h: regular and sparse matrix ///////////////////////////
#ifndef __matrixF_h
#define __matrixF_h
///////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <map>
#include <functional>
#include <iostream>
#include "util.h"
#include "vectorF.h"
using std::vector; using std::less; using std::map; using std::unary_function;
///////////////////////////////////////////////////////////////////////////////////
typedef vector<double> dVec;                   // normal double precision vector //
typedef vector<dVec> dMat;                     // normal double precision matrix //
typedef map<int,double,less<int> > dSpVec;     // sparse double precision vector //
typedef map<int,dSpVec,less<int> > dSpMat;     // sparse double precision matrix //
///////////////////////////////////////////////////////////////////////////////////
class Matrix {
 protected:
  dMat m_;
  dVec column_(int j) const;
 public:
  Matrix():m_() {}
  Matrix(const Matrix& o):m_(o.m_) {}
  Matrix(const dMat& m);
  Matrix(double v,int m,int n);
  int M() const { return m_.size(); }
  int N() const { if(!m_.size()) return 0; else return m_[0].size(); }
  int Rows() const { return m_.size(); } 
  int Cols() const { if(!m_.size()) return 0; else return m_[0].size(); }
  Matrix Product(const Matrix& o) const;
  dVec   Product(const dVec& v) const;
  void Scale(double b);
  void Add(double a);
  dVec operator[](int i) const { return dVec(m_[i]); }
  int size() const { return m_.size(); }
  Matrix Transpose() const;
};
///////////////////////////////////////////////////////////////////////////////////
class SparseMatrix{
 protected:
  dSpMat sm_; int m_,n_;
  dVec column_(int j) const;
  dVec row_(int i) const;
 public:
  SparseMatrix():sm_(),m_(),n_() {}
  SparseMatrix(const SparseMatrix& o):sm_(o.sm_),m_(o.m_),n_(o.n_) {}
  SparseMatrix(const dSpMat& sm,int m,int n):sm_(sm),m_(m),n_(n) {}
  SparseMatrix(const dMat& m);
  SparseMatrix(const Matrix& m);
  int M() const { return m_; }
  int N() const { return n_; }
  int Rows() const { return m_; }
  int Cols() const { return n_; }
  //  SparseMatrix Product(const SparseMatrix& o) const;
  Matrix       Product(const Matrix& o) const;
  dVec         Product(const dVec& v) const;
  void Scale(double b);
  Matrix Transpose() const;
  dVec operator[](int i) const { return dVec(row_(i)); }
};
///////////////////////////////////////////////////////////////////////////////////
class OuterProduct:public binary_function<dVec,dVec,Matrix> {
 public:
  Matrix operator()(const dVec& v1,const dVec& v2) const;
};
///////////////////////////////////////////////////////////////////////////////////
class Identity:public unary_function<int,Matrix> {
 public:
  Matrix operator()(int N) const;
};
///////////////////////////////////////////////////////////////////////////////////
class ZeroesMatrix:public binary_function<int,int,Matrix> {
 public:
  Matrix operator()(int M,int N) const;
};
///////////////////////////////////////////////////////////////////////////////////
class OnesMatrix:public binary_function<int,int,Matrix> {
 public:
  Matrix operator()(int M,int N) const;
};
///////////////////////////////////////////////////////////////////////////////////
class RandomMatrix:public binary_function<int,int,Matrix> {
 public:
  Matrix operator()(int M,int N) const;
};
///////////////////////////////////////////////////////////////////////////////////
dVec operator * (const dVec& v, const Matrix& m);
dVec operator * (const Matrix& m, const dVec& v);
Matrix operator * (const Matrix& m1, const Matrix& m2);
///////////////////////////////////////////////////////////////////////////////////
std::ostream& operator << (std::ostream& os,const Matrix& m);
std::ostream& operator << (std::ostream& os,const SparseMatrix& m);
///////////////////////////////////////////////////////////////////////////////////
#endif
