/////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef __graph_c
#define __graph_c
#include <vector>
typedef std::vector<int> iVec;
typedef std::vector<iVec> iMat;
/////////////////////////////////////////////////////////////////////////////////////////////////
struct Edge {
  int v,w;
  Edge(int V= -1,int W= -1):v(V),w(W) {}
};
/////////////////////////////////////////////////////////////////////////////////////////////////
//class DenseGraph:public GRAPH {
class DenseGraph {
 private:
  int Vcnt,Ecnt; bool digraph; iMat adj;
 public:
  //  DenseGraph():GRAPH() {}
  ~DenseGraph() {}
  DenseGraph(int V,bool Digraph=false):Vcnt(V),Ecnt(0),digraph(Digraph),adj(V) {
    for(int i=0;i<V;++i) adj[i].assign(V,0);
  }
  DenseGraph():Vcnt(),Ecnt(),digraph(),adj() {}
  DenseGraph(const DenseGraph& o):Vcnt(o.Vcnt),Ecnt(o.Ecnt),digraph(o.digraph),adj(o.adj) {}
  int insert(Edge e) {
    int v=e.v, w=e.w;
    if(v>=Vcnt || w>=Vcnt) return Ecnt;
    if(! adj[v][w]) ++Ecnt;
    adj[v][w]= 1;
    if(!digraph) adj[w][v]=1;
    return Ecnt;
  }
  int remove(Edge e) {
    int v=e.v, w=e.w;
    if(v>=Vcnt || w>=Vcnt) return Ecnt;
    if(adj[v][w]) --Ecnt;
    adj[v][w]= 0;
    if(!digraph) adj[w][v]= false;
    return Ecnt;
  } 
  bool edge(int v,int w) const {
    if(v>=Vcnt || w>=Vcnt) return false;
    return (adj[v][w]!=0);
  }
  int V() const { return Vcnt; }
  int E() const { return Ecnt; }
  bool directed() const { return digraph; }
  int vEdges(int v) const { 
    if(v >= Vcnt) return -1; 
    int cnt=0;
    for(int i=0;i<Vcnt;++i) if(adj[v][i]) ++cnt;
    return cnt;
  }
  iVec vEdgeList(int v) const {
    if(v >= Vcnt) return iVec();
    iVec el; el.reserve(Vcnt);
    for(int i=0;i<Vcnt;++i) if(adj[v][i]) el.push_back(i);
    return el;
  }
  int unLinkV(int v) {
    if(v >= Vcnt) return -1;
    for(int i=0;i<Vcnt;++i)
      if(adj[v][i] || adj[i][v]) { 
	--Ecnt; if(digraph) --Ecnt; 
	adj[v][i]=adj[i][v]=0; 
      }
    return Ecnt;
  }
  class adjIterator;
  friend class adjIterator;
};
/////////////////////////////////////////////////////////////////////////////////////////////////
class DenseGraph::adjIterator {
 private:
  const DenseGraph& g;
  int i,v;
 public:
  adjIterator(const DenseGraph& G,int V):g(G),v(V),i(-1) {}
  int beg() {   i= -1;   return nxt(); }
  int nxt() {  
    for(i++; i< g.V(); ++i) if(g.adj[v][i]) return i; 
    return -1;
  }
  bool end() { return i>=g.V(); }
};
/////////////////////////////////////////////////////////////////////////////////////////////////
template <class Graph> class CC {
  const Graph& G;
  int ccnt;
  iVec id;
  void ccR(int w) {
    id[w]=ccnt;
    typename Graph::adjIterator A(G,w);
    for(int v=A.beg(); !A.end(); v=A.nxt())
      if(id[v] == -1) ccR(v);
  }
 public:
  CC(const Graph& g):G(g),ccnt(0),id(G.V(),-1) {
    for(int v=0;v<G.V(); ++v)
      if(id[v] == -1) { ccR(v); ++ccnt; }
  }
  int count() const { return ccnt; }
  bool connect(int s,int t) const {
    if(s>=G.V() || t>=G.V()) return false;
    return id[s]==it[t];
  }
  iVec component(int s) const {
    iVec ret; 
    if(s>=ccnt) return ret;
    int v=G.V();
    ret.reserve(v);
    for(int i=0;i<v;++i) 
      if(id[i]==s) ret.push_back(i);
    return ret;
  }
};
/////////////////////////////////////////////////////////////////////////////////////////////////
template <class Graph> class SEARCH {
 protected:
  const Graph& G;
  int cnt;
  iVec ord;
  virtual void searchC(Edge)=0;
  void search() {
    for (int v=0; v<G.V(); ++v) if (ord[v] == -1) searchC(Edge(v,v)); 
  }
 public:
  SEARCH(const Graph& G):G(G), ord(G.V(), -1), cnt(0) {}
  int operator[](int v) const { return ord[v]; }
};
/////////////////////////////////////////////////////////////////////////////////////////////////
template <class Graph> class EC: public SEARCH<Graph> {
  int bcnt;
  iVec low;
  iMat bridges;
  void searchC(Edge e) {
    int w=e.w;
    ord[w]= cnt++; low[w]= ord[w];
    typename Graph::adjIterator A(G,w);
    for(int t=A.beg(); !A.end(); t=A.nxt())
      if(ord[t] == -1) {
	searchC(Edge(w,t));
	if(low[w] > low[t]) low[w]=low[t];
	if(low[t] == ord[t]) {
	  ++bcnt; //w-t is a bridge
	  iVec tmp(2,0);
	  if(w<t) { tmp[0]=w; tmp[1]=t; }
	  else    { tmp[1]=w; tmp[0]=t; }
	  bridges.push_back(tmp);
	}
      } else if (t != e.v) 
	if(low[w] > ord[t]) low[w]= ord[t];
  }
 public:
  EC(const Graph& G): SEARCH<Graph>(G),bcnt(0),low(G.V(),-1) { search(); }
  int count() const { return bcnt; }
  iVec Bridge(int i) const { if(i<bcnt) return bridges[i]; return iVec(); }
  iMat Bridges() const { return bridges; }
};
/////////////////////////////////////////////////////////////////////////////////////////////////
#endif
/////////////////////////////////////////////////////////////////////////////////////////////////
