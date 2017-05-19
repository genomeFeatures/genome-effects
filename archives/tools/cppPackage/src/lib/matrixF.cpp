///////////////////////////////////////////////////////////////////////////////////
#include "../lib/logFile.h"
#include "../lib/matrixF.h"
#include <time.h>
///////////////////////////////////////////////////////////////////////////////////
#include <string>
using std::string;
///////////////////////////////////////////////////////////////////////////////////
Matrix::Matrix(double v,int m,int n):m_() {
  if(m<=0) { WriteIntToErr()(string("invalid m in Matrix constructor"),m); return; }
  if(n<=0) { WriteIntToErr()(string("invalid n in Matrix constructor"),n); return; }
  dVec dummy; dummy.reserve(n);
  for(int i=0;i<n;i++) dummy.push_back(v);
  m_.reserve(m);
  for(int i=0;i<m;i++) m_.push_back(dummy);
  return;
}
///////////////////////////////////////////////////////////////////////////////////
Matrix::Matrix(const dMat& m):m_(m) {
  if(!m_.size()) { WriteToErr()(string("zero rows in dMat m in Matrix constructor")); return; }
  if(!m_[0].size()) { WriteToErr()(string("zero cols in dMat m in Matrix constructor")); return; }
  for(int i=1;i<m_.size();i++)
    if(m_[i].size()!=m_[0].size()) {
      WriteToErr()(string("rows have varying sizes in dMat m in Matrix constructor")); return; }
  return;
}
///////////////////////////////////////////////////////////////////////////////////
dVec Matrix::column_(int j) const {
  if(j>=m_[0].size()) {
    WriteIntToErr()(string("index j too high in Matrix::column(int j)"),j);
    return dVec();
  }
  dVec c; c.reserve(m_.size());
  for(int i=0;i<m_.size();i++) c.push_back(m_[i][j]);
  return dVec(c);
}
///////////////////////////////////////////////////////////////////////////////////
Matrix Matrix::Product(const Matrix& o) const {
  if(m_[0].size()!=o.M()) {
    WriteToErr()(string("Matrix1 cols != Matrix2 rows in Matrix::Product"));
    return Matrix();
  }
  int m=m_.size(),n=m_[0].size(),n2=o.N();
  dMat p; p.reserve(m);
  for(int i=0;i<m;i++) {
    dVec row_; row_.reserve(n2);
    for(int i2=0;i2<n2;i2++) {
      double tsum=0.0;
      for(int j=0;j<n;j++) tsum += m_[i][j]*o.m_[j][i2];
      row_.push_back(tsum);
    }
    p.push_back(row_);
  }
  return Matrix(p);
}
///////////////////////////////////////////////////////////////////////////////////
dVec   Matrix::Product(const dVec& v) const {
  if(m_[0].size()!=v.size()) {
    WriteToErr()(string("Matrix cols != vector rows in Matrix::Product"));
    return dVec();
  }
  int m=m_.size(),n=m_[0].size();
  dVec p; p.reserve(m);
  for(int i=0;i<m;i++) {
    double tsum=0.0;
    for(int j=0;j<n;j++) tsum += m_[i][j]*v[j];
    p.push_back(tsum);
  }
  return dVec(p);
}
///////////////////////////////////////////////////////////////////////////////////
void Matrix::Scale(double b) {
  int m=m_.size(),n=m_[0].size();
  for(int i=0;i<m;i++)
    for(int j=0;j<n;j++)
      m_[i][j] *= b;
  return;
}
///////////////////////////////////////////////////////////////////////////////////
void Matrix::Add(double a) {
  int m=m_.size(),n=m_[0].size();
  for(int i=0;i<m;i++)
    for(int j=0;j<n;j++)
      m_[i][j] += a;
  return;
}
///////////////////////////////////////////////////////////////////////////////////
Matrix Matrix::Transpose() const {
  dMat m;
  for(unsigned i=0;i<m_[0].size();i++)
    m.push_back(column_(i));
  return Matrix(m);
}
///////////////////////////////////////////////////////////////////////////////////
SparseMatrix::SparseMatrix(const dMat& m):sm_(),m_(),n_() {
  if(!m.size()) {
    WriteToErr()(string("zero rows in dMat m in SparseMatrix constructor"));
    return;
  }
  if(!m[0].size()) {
    WriteToErr()(string("zero columns in dMat m in SparseMatrix constructor"));
    return;
  }
  int M=m.size(),N=m[0].size();
  for(int i=1;i<M;i++)
    if(m[i].size()!=N) {
      WriteToErr()(string("diff sized rows in dMat m in SparseMatrix constructor"));
      return;
    }
  m_=M; n_=N;
  for(int i=0;i<m_;i++) {
    dSpVec row;
    for(int j=0;j<N;j++) if(m[i][j]!=0.0) row[j]= m[i][j];
    if(row.size()) sm_[i]= row;
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////
SparseMatrix::SparseMatrix(const Matrix& m):sm_(),m_(m.M()),n_(m.N()) {
  if(!m.M()) {
    WriteToErr()(string("zero rows in Matrix m in SparseMatrix constructor"));
    return; 
  }
  if(!m.N()) {
    WriteToErr()(string("zero columns in Matrix m in SparseMatrix constructor"));
    return;
  }
  for(unsigned i=0;i<m_;i++) {
    dSpVec row;
    for(int j=0;j<n_;j++) if(m[i][j]!=0.0) row[j]= m[i][j];
    if(row.size()) sm_[i]= row;
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////
dVec SparseMatrix::column_(int j) const {
  if((j<0)||(j>=n_)) {
    WriteIntToErr()(string("invalid j in SparseMatrix::column_: "),j);
    return dVec();
  }
  dVec v;
  for(int i=0;i<m_;i++) {
    dSpMat::const_iterator r= sm_.find(i);
    if(r==sm_.end())                           v.push_back(0.0); 
    else {
      dSpVec::const_iterator c= (*r).second.find(j);
      if(c==(*r).second.end())                 v.push_back(0.0);
      else                                     v.push_back((*c).second);
    }
  }
  return dVec(v);
}
///////////////////////////////////////////////////////////////////////////////////
dVec SparseMatrix::row_(int i) const {
  if((i<0)||(i>=m_)) {
    WriteIntToErr()(string("invalid i in SparseMatrix::column_: "),i);
    return dVec();
  }
  dVec v;
  dSpMat::const_iterator r= sm_.find(i);
  for(int j=0;j<n_;j++)
    if(r==sm_.end())                      v.push_back(0.0);
    else {
      dSpVec::const_iterator c= (*r).second.find(j);
      if(c==(*r).second.end())            v.push_back(0.0);
      else                                v.push_back((*c).second);
    }
  return dVec(v);
}
///////////////////////////////////////////////////////////////////////////////////
//SparseMatrix SparseMatrix::Product(const SparseMatrix& o) const {
//  return SparseMatrix();
//}
///////////////////////////////////////////////////////////////////////////////////
Matrix SparseMatrix::Product(const Matrix& o) const {
  if(n_!=o.M()) {
    WriteToErr()(string("Matrix and SparseMatrix sizes incompatible in SparseMatrix product"));
    return Matrix();
  }
  dMat p; p.reserve(m_);
  for(dSpMat::const_iterator i=sm_.begin();i!=sm_.end();i++) {
    dVec row; row.reserve(o.N());
    for(unsigned i2=0;i2<o.N();i2++) {
      double tsum=0.0;
      for(dSpVec::const_iterator j=(*i).second.begin();j!=(*i).second.end();j++)
	tsum += (*j).second * o[(*i).first][(*j).first];
      row.push_back(tsum);
    }
    p.push_back(row);
  }
  return Matrix(p);
}
///////////////////////////////////////////////////////////////////////////////////
dVec SparseMatrix::Product(const dVec& v) const {
  if(n_!=v.size()) {
    WriteToErr()(string("dVec and SparseMatrix sizes incompatible in SparseMatrix product"));
    return dVec();
  }
  dVec p; p.reserve(m_);
  for(dSpMat::const_iterator i=sm_.begin();i!=sm_.end();i++) {
    double tsum=0.0;
    for(dSpVec::const_iterator j=(*i).second.begin();j!=(*i).second.end();j++)
      tsum += (*j).second * v[(*j).first];
    p.push_back(tsum);
  }
  return dVec(p);
}
///////////////////////////////////////////////////////////////////////////////////
void SparseMatrix::Scale(double b) {
  for(dSpMat::iterator i=sm_.begin();i!=sm_.end();i++)
    for(dSpVec::iterator j= (*i).second.begin();j!=(*i).second.end();j++)
      (*j).second *= b;
}
///////////////////////////////////////////////////////////////////////////////////
Matrix SparseMatrix::Transpose() const {
  dMat m;
  for(unsigned i=0;i<n_;i++)
    m.push_back(column_(i));
  return Matrix(m);
}
///////////////////////////////////////////////////////////////////////////////////
Matrix OuterProduct::operator()(const dVec& v1,const dVec& v2) const {
  dMat m;
  for(int i=0;i<v1.size();i++) {
    dVec row;
    for(int j=0;j<v2.size();j++) row.push_back(v1[i]*v2[j]);
    m.push_back(row);
  }
  return(Matrix(m));
}
///////////////////////////////////////////////////////////////////////////////////
Matrix Identity::operator()(int N) const {
  dMat m;
  for(unsigned i=0;i<N;i++) {
    dVec row;
    for(unsigned j=0;j<N;j++)
      if(i==j) row.push_back(1.0);
      else     row.push_back(0.0);
    m.push_back(row);
  }
  return Matrix(m);
}
///////////////////////////////////////////////////////////////////////////////////
Matrix OnesMatrix::operator()(int M,int N) const {
  dMat m;
  for(unsigned i=0;i<M;i++) {
    dVec row;
    for(unsigned j=0;j<N;j++) row.push_back(1.0);
    m.push_back(row);
  }
  return Matrix(m);
}
///////////////////////////////////////////////////////////////////////////////////
Matrix ZeroesMatrix::operator()(int M,int N) const {
  dMat m;
  for(unsigned i=0;i<M;i++) {
    dVec row;
    for(unsigned j=0;j<N;j++) row.push_back(0.0);
    m.push_back(row);
  }
  return Matrix(m);
}
///////////////////////////////////////////////////////////////////////////////////
Matrix RandomMatrix::operator()(int M,int N) const {
  dMat m; 
  long seed; time_t t0; time(&t0); seed= -1*t0;
  for(unsigned i=0;i<M;i++) {
    dVec r;
    for(unsigned i=0;i<N;i++) r.push_back(rand3()(&seed));
    m.push_back(r);
  }
  return Matrix(m);
}
///////////////////////////////////////////////////////////////////////////////////
dVec operator * (const dVec& v, const Matrix& m) {
  if(m.M()!=v.size()) {
    WriteToErr()(string("Matrix rows != vector elements in operator * (dVec, Matrix)"));
    return dVec();
  }
  int M=m.N(),N=m.M();
  dVec p; p.reserve(M);
  for(int i=0;i<M;i++) {
    double tsum=0.0;
    for(int j=0;j<N;j++) tsum += m[j][i]*v[j];
    p.push_back(tsum);
  }
  return dVec(p);
}
///////////////////////////////////////////////////////////////////////////////////
dVec operator * (const Matrix& m, const dVec& v) {
  if(m.N()!=v.size()) {
    WriteToErr()(string("Matrix cols != vector elements in operator * (Matrix, dVec)"));
    return dVec();
  }
  int M=m.M(),N=m.N();
  dVec p; p.reserve(M);
  for(int i=0;i<M;i++) {
    double tsum=0.0;
    for(int j=0;j<N;j++) tsum += m[i][j]*v[j];
    p.push_back(tsum);
  }
  return dVec(p);
}
///////////////////////////////////////////////////////////////////////////////////
Matrix operator * (const Matrix& m1,const Matrix& m2) {
  if(m1.N()!=m2.M()) {
    WriteToErr()(string("Matrix1 cols != Matrix2 rows in Matrix::Product"));
    return Matrix();
  }
  int m=m1.M(),n=m1.N(),n2=m2.N();
  dMat p; p.reserve(m);
  for(int i=0;i<m;i++) {
    dVec row_; row_.reserve(n2);
    for(int i2=0;i2<n2;i2++) {
      double tsum=0.0;
      for(int j=0;j<n;j++) tsum += m1[i][j]*m2[j][i2];
      row_.push_back(tsum);
    }
    p.push_back(row_);
  }
  return Matrix(p);
}
///////////////////////////////////////////////////////////////////////////////////
std::ostream& operator << (std::ostream& os,const Matrix& m) {
  for(unsigned i=0;i<m.M();i++) os << m[i];
  return os;
}
///////////////////////////////////////////////////////////////////////////////////
std::ostream& operator << (std::ostream& os,const SparseMatrix& m) {
  for(unsigned i=0;i<m.M();i++) os << m[i];
  return os;
}
///////////////////////////////////////////////////////////////////////////////////
