/*  Copyright (C) 2000, 2006  Marc Toussaint (mtoussai@inf.ed.ac.uk)
    under the terms of the GNU LGPL (http://www.gnu.org/copyleft/lesser.html)
    see the `std.h' file for a full copyright statement  */

#ifndef MT_list_h
#define MT_list_h

//#define MT_LISTCHECK

#include "std.h"
#include "array.h"
#include <algorithm>
#define INSbefore 1
#define INSafter 0

#ifndef forNodes
#  define forNodes(n,L)            for(n=L.first;n;n=n->next)
#  define forNodes_rev(n,L)        for(n=L.last;n;n=n->prev)
#  define forNodes_save(n,m,L)     for(n=L.first,m=n?n->next:0; n; n=m,m=n?n->next:0)
#endif
#ifndef for_nodes
#  define for_nodes(n) for(n=first;n;n=n->next)
#  define for_nodes_rev(n) for(n=last;n;n=n->prev)
#  define for_nodes_save(n,sn) for(n=first,sn=n?n->next:0; n; n=sn,sn=n?n->next:0)
//#  define for_nodes_rev(n,sn) for(n=last;n;n=n->prev)
#endif

//===========================================================================
//
// List class
//

namespace MT{
/*! A general list object.*/
template<class T>
class List{
public: //types
  typedef List<T> ListT;

  class nodeClass:public T{
  public:
    nodeClass *next,*prev;
    uint index; //gives the node's index in the list's \c nodes array
    double w;
    nodeClass():T(){ setNull(); }
    nodeClass(const T& n):T(n){ setNull(); }
    T& data(){ return (T&)(*this); }
  private:
    void setNull(){ next=prev=0; w=0.; }
  };
  /*! A node is a pointer to a \c nodeClass object which is derived
      from your data type \c T; see List class description. */ 
  typedef nodeClass *node; //the iterator and pointer

  /*! this class type realizes a mapping from each node of the list
      to a datum of type \c T */
  template<class D>
  class node_map{
  public:
    MT::Array<D> A;
    node_map(){}
    node_map(const ListT& L){ init(L); }
    void init(const ListT& L){ A.resizeCopy(L.N); }
    D& operator[](node n){ return A(n->index); }
  };

  typedef bool (*NodeCompare)(const T& a,const T& b);

public:
  uint N; //!< number of nodes
  /*! the first node of the list */
  node first;
  /*! the last node of the list */
  node last;

  /*! an array of all nodes (pointers to nodeClasses): Actually all
    nodes are linked in a list via their prev and next pointers -- in
    principle this would be enough to keep track of them. Still, this
    array additionally comprises pointers to all nodes. (Resize is
    quick because of the MT::Array::flexiMem and MT::Array::memMove
    flags turned on.) This allows faster random access and quick
    checking if a node belongs to this list. Note though that the
    ordering defined by the linked list of nodes, and the indexing
    defined by this array, generally don't coincide. Use
    indexAsOrdered() or orderAsIndexed() to synchronize the ordering
    and indexing.*/
  MT::Array<node> nodes;

  /*! tells whether the list was defragmented (see defragment()) */
  bool defragmented;

public:
  List(){
    N=0;
    first=last=0;
    nodes.memMove=true;
    defragmented=false;
  }
  List(const List& L){
    N=0;
    first=last=0;
    nodes.memMove=true;
    defragmented=false;
    copy(L);
  }
  ~List(){
    clear();
  }


public://{ basic access:
  /*! the i-th node */
  node operator()(uint i){ return nodes(i); }
  /*! a randomly (uniformly) chosen node */
  node rndNode(){ if(!N) return 0; return nodes(rnd(N)); }
  /*! return true if n is a node of this list */
  bool checkNode(node n) const{ return (n && n->index<N && nodes(n->index)==n); }


public://{ the core methods: new and delete nodes
  /*! inserts a new node after or before \c where (default: at the end
      of the node list); opt has to be INSafter or INSbefore,
      respectively */
  node new_node(node where=0,int opt=INSbefore,const T* x=0){
    CHECK(opt==INSafter || opt==INSbefore,"List::new_node");
    CHECK(!defragmented,"new_node: can't mutate defragmented list!");
#ifdef MT_LISTCHECK
    checkValid();
#endif
    if(!where){
      if(opt==INSbefore) where=last; else where=0;
    }else{
      CHECK(checkNode(where),"new_node: where is not a node");
      if(opt==INSbefore) where=where->prev;
    }
    node n;
    if(!x) n=new nodeClass;
    else   n=new nodeClass(*x);
    if(!n) HALT("MT::List failed node memory allocation (list size="<<N*sizeof(nodeClass)<<"bytes)");
    //link up in list:
    if(!N){
      n->prev=0;      n->next=0;
      first=n;	      last=n;
    }else{
      if(where){ n->prev=where; n->next=where->next; }else{ n->prev=0; n->next=first; }
      if(n->prev) n->prev->next=n; else first=n;
      if(n->next) n->next->prev=n; else last=n;
    }
    //register in table:
    n->index=N;
    N++;
    nodes.resizeCopy(nodes.N+1);
    nodes(n->index)=n;

#ifdef MT_LISTCHECK
    checkValid();
#endif
    return n;
  }
  /*! appends a new node and initialize it by \c x */
  node new_node(const T& x){ return new_node(last,INSafter,&x); }
  /*! pushs x in front of the list */
  node add(const T& x){ return new_node(first,INSbefore,&x); }
  /*! appends x to the list */
  node append(const T& x){ return new_node(last,INSafter,&x); }

  /*! deletes the node. \c n becomes invalid! */
  void del_node(node n){
    CHECK(checkNode(n),"del_node: is not a node");
    CHECK(!defragmented,"del_node: can't mutate defragmented list!");
#ifdef MT_LISTCHECK
    checkValid();
#endif
    if(!n) return;
    //unlink from list:
    if(n==first) first=n->next; else n->prev->next=n->next;
    if(n==last)  last =n->prev; else n->next->prev=n->prev;
    //delete from table:
    nodes(n->index)=nodes(N-1);
    nodes(n->index)->index=n->index;
    N--;
    nodes.resizeCopy(nodes.N-1);
    delete n;
#ifdef MT_LISTCHECK
    checkValid();
#endif
  }

  //! same as del_node(n)
  void remove(const T& x){
    node n,ns;
    for_nodes_save(n,ns) if(n->data()==x){ del_node(n); break; }
  }

  /*! deletes all nodes */
  void clear(){
    uint i;
    if(!defragmented){
      for(i=0;i<N;i++) delete nodes(i);
    }else{
      if(N)  delete[] nodes(0);
    }
    first=last=0;
    N=0;
    nodes.resize(0);
#ifdef MT_LISTCHECK
    checkValid();
#endif
  }


public://{ standard looping

  /*! for(First(n);n;Next(n)) instead of for(n=first;n;n=n->next) */
  void First(node& n) const{ n=first; }
  /*! .. */
  void Next(node& n) const{ n=n?n->next:first; }
  bool CheckNext(node& n) const{ return n!=0; }
  /*! for(Last(n);n;Prev(n)) instead of for(n=last;n;n=n->prev) */
  void Last(node& n) const{ n=last; }
  /*! .. */
  void Prev(node& n) const{ n=n?n->prev:last; }
  bool CheckPrev(node& n) const{ return n!=0; }


public://{ list assignment
  /*! copy operator */
  ListT& operator=(const ListT& L){
    copy(L);
    return *this;
  }

  /*! *this becomes a copy of L */
  void copy(const ListT& L){
    node v;
    clear();
    forall(v,L) new_node(*v);
  }

  void set(uint n,const T& x){
    clear();
    for(uint i=0;i<n;i++) new_node(x);
  }

  /*! inserts L in *this before or after the node l;
      copies elements instead of stealing (L remains untouched) */
  void insert(const ListT& L,node where=0,int opt=INSbefore){
    node v;
    if(opt==INSbefore) forall(v,L) new_node(where,INSbefore,v);
    else           forall_rev(v,L) new_node(where,INSafter,v);
  }

  /*! steals all elements from L; L becomes empty!; the indexing of *this is reset */
  void unite(ListT& L,node where=0,int opt=INSbefore){
    CHECK(!defragmented && !L.defragmented,"can't join defragmented lists");
    if(!L.N) return;
    node before,after;
    if(opt==INSbefore){
      if(where){ before=where->prev; after=where; }
      else{ before=last; after=0; }
    }else{
      if(where){ before=where; after=where->next; }
      else{ before=0; after=first; }
    }

    if(before){ before->next=L.first; L.first->prev=before; }else first=L.first;
    if(after){  after->prev =L.last;  L.last->next =after;  }else last=L.last;

    N=N+L.N;
    nodes.resize(nodes.N+L.N);
    uint i=0; node n;
    for_nodes(n){ n->index=i; nodes(i)=n; i++; }
    CHECK(i==N,"dammit");

    L.N=0;
    L.nodes.resize(0);
    L.first=L.last=0;

#ifdef MT_LISTCHECK
    checkValid();
    L.checkValid();
#endif
  }

  /*! copy the list L */
  template<class FList,class CopyStruct>
  void setForeignList(const FList& L,const CopyStruct&){ 
    node nn;
    typename FList::node n;
    clear();
    forall(n,L){
      nn=new_node();
      CopyStruct::copy(*nn,*n);
    }
  }

  /*! generates a list of \c n new nodes */
  void setEmptySize(uint n){
    CHECK(n>0,"");
    clear();
    for(uint i=0;i<n;i++) new_node();
  }

  void resize(uint n){
    if(N>n){ while(N>n) del_node(last); }
    if(N<n){ while(N<n) new_node(); }
  }

public://{ list comparison
  bool operator==(const ListT& L){
    if(N!=L.N) return false;
    node n,m;
    for(n=first,m=L.first;n && m;n=n->next,m=m->next){
      if(!(n->data()==m->data())) return false;
    }
    return true;
  }

public://{ list manipulation
  /*! duplicates n */
  node duplicate(node n){
    CHECK(checkNode(n),"List::duplicate: invalid node");
    node m;
    m=new_node(n,INSafter,n);
    return m;
  }
  //}

private:
  void allocateDefragMem(uint nn){
    clear();
    if(!nn) return;
    uint i;
    node NA;
    NA=new nodeClass [nn];
    for(i=0;i<nn;i++){
      NA[i].index=i;
      if(i+1<nn) NA[i].next=&NA[i+1]; else NA[i].next=0;
      if(i>0)    NA[i].prev=&NA[i-1]; else NA[i].prev=0;
    }
    first=&NA[0];
    last =&NA[nn-1];
    N=nn;
    nodes.resize(nn);

    for(i=0;i<nn;i++) nodes(i)=&NA[i];

    defragmented=true;
  }



public://{ ordering and indexing
  /*! \c this becomes a defragmented copy of the list \c L. This means
      that all the memory for nodes is allocated in such a
      way that they become linearly aligned in ordinary C arrays. This
      alignment is w.r.t.\ to the INDEXING of the list, not the
      ordering. Actually, you can do your own indexing as you want to
      specify the alignment, or use List::indexAsOrdered.

      In particular on Intel architectures, this makes recursion much
      faster (Up to 5 times faster for large lists! Maybe because the
      memory can be uploaded to the cache in huge
      blocks?). E.g. recursion over all nodes can now be done by
      ordinary pointer increments. Note that this is a one-way ticket:
      node deletion/insertion are not possible for a defragmented
      list and you cannot un-defragment a list directly. You can
      though copy it into another list, which is then un-defragmented
      (and then, if you like, copy it back into this list) */
  void defragment(const List& L){
#ifdef MT_LISTCHECK
    checkValid();
#endif
    uint i,nn=L.N;
    node n,Fi,La;
    nodeClass *NA=0;
    if(nn) NA=new nodeClass [nn];
    for(i=0;i<nn;i++){
      NA[i]=*L.nodes(i);
      if(L.nodes(i)->next)     NA[i].next    =&NA[L.nodes(i)->next->index];
      if(L.nodes(i)->prev)     NA[i].prev    =&NA[L.nodes(i)->prev->index];
    }

    if(L.first) Fi=&NA[L.first->index]; else Fi=0;
    if(L.last)  La=&NA[L.last->index];  else La=0;

    clear();

    first=Fi;
    last=La;
    N=nn;
    nodes.resize(nn);
    for_nodes(n){
      nodes(n->index)=n;
    }
    defragmented=true;

#ifdef MT_LISTCHECK
    checkValid();
#endif
  }

  /*! sorts the indexing (i.e. the List::nodes array) according to
      the ordering (as given by the node list); afterward, the
      ordering and the indexing coincide */
  void indexAsOrdered(){
    node n; uint i=0;
    for_nodes(n){ n->index=i; nodes(i)=n; i++; }
#ifdef MT_LISTCHECK
    checkValid();
#endif
  }

  //! permutes the indices of the i-th and j-th node
  void permuteIndex(uint i,uint j){
    node n=nodes(i); nodes(i)=nodes(j); nodes(j)=n;
    nodes(i)->index=i; nodes(j)->index=j;
#ifdef MT_LISTCHECK
    checkValid();
#endif
  }

  //! indexes all nodes in random order
  void randomIndexing(){
    MT::Permutation pi;
    pi.random(N);
    node n;
    uint i=0;
    for_nodes(n){ n->index=i; nodes(i)=n; i++; }
#ifdef MT_LISTCHECK
    checkValid();
#endif
  }

  /*! this function reliks the node list (i.e., changes the actual
      ordering of the nodes!) according to the order of the node
      indices (i.e., the List::nodes array); e.g., after sortIndexing() or
      if you permuted some indices */
  void orderAsIndexed(){
    uint i;
    if(N){
      first=nodes(0);    first->prev=0; 
      last=nodes(N-1); last->next=0; 
      if(N>1){ first->next=nodes(1); last->prev=nodes(N-2); }
      else{ first->next=0; last->prev=0; }
      for(i=1;i<N-1;i++){
	nodes(i)->prev=nodes(i-1);
	nodes(i)->next=nodes(i+1);
      }
    }
#ifdef MT_LISTCHECK
    checkValid();
#endif
  }

  /*! sorts the node indexing according to postive comparison operator */
  template<class CompareOperator>
  void sortIndexing(CompareOperator comp){
    uint i;
    node *a=(node*)nodes.p,*b=(node*)nodes.pstop;
    std::sort(a,b,comp);
    for(i=0;i<nodes.N;i++) nodes(i)->index=i;
  }

  struct LowerOperator{
    bool operator()(const node& a,const node& b) const{ return a->data() < b->data(); }
  };
  //! same as sortIndexingOp(LowerOperator()); 
  void sortIndexing(){ sortIndexingOp(LowerOperator()); }

  //! same as sortIndexing(); orderAsIndexed();
  template<class CompareOperator>
  void sort(CompareOperator comp){ sortIndexing(comp); orderAsIndexed(); }

  void sort(){ sortIndexing(LowerOperator()); orderAsIndexed(); }

  /*! inserts x in a sorted list (from small to large)
      so that it remains sorted, using findInSorted.
      WARNING: this breaks the data structure,
      namely the indices do not match anymore the indexing
      (n!=nodex(n->index)) */
  node addSorted(const T& x,NodeCompare comp){
    if(!N) return new_node(0,INSafter,&x);
    uint i=findInSorted(x,comp);
    node n;
    if(i==0)
      n=new_node(0,INSafter,&x);
    else{
      nodes(i-1)->index=i-1;
      n=new_node(nodes(i-1),INSafter,&x);
    }
    if(nodes.N>i+1) memmove(&nodes(i+1),&nodes(i),sizeof(uint)*(nodes.N-i-1));
    nodes(i) = n;
    CHECK(nodes(0)==first && nodes(N-1)==last,"hack");
    return n;
  }
  void delSorted(node n,NodeCompare comp){
    if(n==last){ last->index=N-1; del_node(last); return; }
    uint i=findInSorted(n->data(),comp);
    for(;;){
      if(n==nodes(i)){ delSorted(i); return; }
      i--;
    }
    HALT("you shouldn't be here!");
  }
  void delSorted(uint i){
    nodes(i)->index=i;
    del_node(nodes(i));
    if(nodes.N>i+1) memmove(&nodes(i),&nodes(i+1),sizeof(uint)*(nodes.N-i-1));
    nodes(N-1)=last;
    CHECK(nodes(0)==first && nodes(N-1)==last,"hack");
  }

  /*! finds the node of a sorted list which is smaller or equal
      than x; takes log N time */
  uint findInSorted(const T& x,NodeCompare comp){
    if(!N) return 0;
    if(comp(x,first->data())) return 0;
    if(!comp(x,last->data())) return N;
    uint lo=0,up=N-1,mi;
    for(;;){
      if(lo+1==up) return up;
      mi=lo+(up-lo)/2;
      if(comp(x,nodes(mi)->data())) up=mi; else lo=mi;
    }
    HALT("you shouldn't be here!");
    return 0;
  }


  /*! returns the ordering of the nodes, i.e., their position in the
    node list */
  void getOrder(node_map<uint>& order) const{
    node n; uint i=0;
    for_nodes(n){ order[n]=i; i++; }
  }

  /*! only for debugging : checks internal lists and tables of nodes
      for consistency */
  bool checkValid() const{
    node n; uint i;

    CHECK(N==nodes.N,"");

    //check node table:
    for(i=0;i<N;i++) CHECK(nodes(i)->index==i,"");

    //check node list:
    uint num=0;
    for_nodes(n){
      if(n->prev==0){ CHECK(first==n,""); }else CHECK(n->prev->next==n,"");
      if(n->next==0){ CHECK(last ==n,""); }else CHECK(n->next->prev==n,"");
      num++;
    }
    CHECK(num==N,"");

    return true;
  }


public://{ some operations - e.g., related to the list being a probabilistic sample set
  /*! resamples the list using Stochastic Universal Sampling; assumes that to each
      element of the list an probability can be associated
      (via the ProbabilityOp) */
  template<class ProbabilityOp>
  void resampleSUS(uint nSamples,const ProbabilityOp& prob){
    double increment=1./nSamples;
    node n,f=first;
    double sum=0,ptr=MT::rnd.uni(0.,increment);
    //sample new elements -- push them in front
    for(n=first;n;n=n->next){
      sum+=prob.P(n);
      while(sum>ptr){ add(n->data()); ptr+=increment; }
    }
    //delete old elements
    f=f->prev; //last new element
    while(f->next) del_node(f->next);
    CHECK(N==nSamples,"resampling of list failed - probabilities not normalized?");
  }

  //! assigns a weight 1/N to each sample
  void weightEqual(){ node n; for_nodes(n) n->w=1./N; }

  //! return the average weighting
  double getAverage(){
    if(!N) return 0.;
    node n;
    double Z=0.;
    for_nodes(n) Z+=n->w;
    return Z/N;
  }

  /*! return the best (highest weighted) sequence and the average
    weighting */
  void getAverageAndBest(double& average,node& best){
    if(!N){ average=0.; return; }
    node n;
    double Z=0.;
    best=0;
    for_nodes(n){ Z+=n->w; if(!best || n->w>best->w) best=n; }
    average=Z/N;
  }

  void checkNorm(){
    node n;
    uint i,mL=0,mS=0;
    double Z=0.;
    for_nodes(n){
      if(n->N>mL) mL=n->N;
      for(i=0;i<n->N;i++){
	CHECK(n->elem(i)>=0,"");
	if(n->elem(i)>(int)mS) mS=n->elem(i);
      }
      Z+=n->w;
    }
    CHECK(fabs(Z-1.)<1e-6,"WeightedList check failed");
    CHECK(mS<100,"WeightedList check failed");
    CHECK(mL<1000,"WeightedList check failed");
  }

  //! normalize the weightings (-> probabilties)
  void normalize(){
    node n;
    double Z=0.;
    for_nodes(n) Z+=n->w;
    if(Z) for_nodes(n) n->w/=Z; else for_nodes(n) n->w=1./N;
  }

  /*! deletes all but the best \c k sequences from the samples
      (weightings remain unchanged (not normalized)) */
  void selectBest(uint k){
    CHECK(k<=N,"can't select best "<<k<<" out of "<<N);
    randomIndexing();
    sortWeights();
    //for(uint i=0;i<k;i++) nodes(i)->w=1./k;
    while(N>k) del_node(nodes(k));
  }

  /*! copies the best \c k sequences from \c P to \c *this */
  void copyBest(ListT& P,uint k){
    if(k>P.N) k=P.N;
    //CHECK(k<=P.N,"can't select best "<<k<<" out of "<<P.N);
    node n;
    clear();
    P.randomIndexing();
    P.sortWeights();
    for(uint i=0;i<k;i++){ n=append(P.nodes(i)->data()); n->w=P.nodes(i)->w; }
  }

  //! weighting become a Boltzmann distributions of themselves
  void Boltzmann(double beta){
    node n;
    for_nodes(n) n->w = exp(beta*n->w);
    normalize();
  }

  struct WeightCmp{ bool operator()(node a,node b){ return a->w > b->w; } };
  struct LexicalCmp{ bool operator()(node a,node b){ return a->data() < b->data(); } };
  struct intA{ static inline double& P(node a){ return a->w; } };

  //! sort list to descending weights (using the STL sort algorithm) 
  void sortWeights(){ sortIndexing(WeightCmp()); }
  //! sort list to increasing lexical order of the sequences
  void sortLexical(){ sortIndexing(LexicalCmp()); }

  /*! Afterwards, the samples will contain \c n samples (generated
      by stochastic universal sampling, or ``roulette wheel
      sampling'') from the distribution it represented before. While
      before weights might have been non-uniform, afterwards all
      samples have equal weight. Effecively, samples of high weight
      are duplicated (or multiplicated) and those of low weight are
      discarded. */
  void resample(uint n){ normalize(); resampleSUS(n,intA()); weightEqual(); }

  //! return the number of different sequences in the sample (which is <=N)
  uint disjoint(){
    uint i,n;
    sortLexical();
    for(i=1,n=1;i<N;i++)
      if(nodes(i)->data()!=nodes(i-1)->data()) n++;
    return n;
  }

  //! deletes multiply occuring sequences (multiset -> set)
  void makeDisjoint(){
    uint i;
    node n,sn;
    sortLexical();
    for(i=1;i<N;i++){
      if(nodes(i)->data()==nodes(i-1)->data()) nodes(i)->w=-1.;
    }
    for_nodes_save(n,sn) if(n->w==-1.) del_node(n);
  }


public://{ I/O
  /*! prototype for \c operator<< */
  void write(std::ostream& os,char* SEP="\n") const{
    node n;
    if(!MT::IOraw) os <<"<List:" <<N <<' ';
    for_nodes(n) os <<SEP <<(*n);
    //for(uint i=0;i<N;i++) os <<SEP <<*nodes(i) <<' ' <<i;
    if(!MT::IOraw) os <<">" <<std::endl;
  }

  /*! prototype for \c operator>> */
  void read(std::istream& is,bool raw=false){
    uint i,nn;
    if(!raw){
      is >>"<List:" >>nn;
      CHECK(is.good(),"could not read list tag");
      clear();
      if(!defragmented) setEmptySize(nn); else allocateDefragMem(nn);
      for(i=0;i<nn;i++){
	is >>*nodes(i);
	CHECK(is.good(),"could not read the "<<i<<"th node of list");
      }
      is >>">";
      CHECK(is.good(),"could not read list stop tag");
    }else{
      clear();
      T x;
      for(;;){
        if(!is.good()) break;
	is >>x;
        if(is.fail()) HALT("error when reading "<<nodes.N<<"th element in readRaw");
        append(x);
      }
    }
  }

  //! opens file to read from 
  void readRaw(const char* filename){
    std::ifstream is; MT::open(is,filename); read(is,true);
  }
  //}
};
}
template<class T>
inline std::istream& operator>>(std::istream& is,MT::List<T>& x){ x.read(is);return is; }
template<class T>
inline std::ostream& operator<<(std::ostream& os,const MT::List<T>& x){ x.write(os); return os; }


//===========================================================================
//
// primitive type list
//

namespace MT{
/*! This class allows to use the List class with primitive types
  (i.e., non-structs/classes like int, double, etc). You have to refer
  to the node data via n->x */
template<class T> class list:public MT::List<basetype<T> >{
public:
  typedef typename MT::List<basetype<T> > ListT;
  typedef typename ListT::node node;
  node new_node(node where=0,int opt=INSafter,const T* x=0){
    if(!x) return ListT::new_node(where,opt,0);
    return ListT::new_node(where,opt,&basetype<T>(*x));
  }
  node new_node(const T& x){ return ListT::new_node(basetype<T>(x)); }
};
}

typedef MT::List<MT::Int> intL;
typedef MT::List<intA> intAL;
typedef MT::List<doubleA> doubleAL;

void random(doubleAL& L,uint N,uint d,double lo,double hi);
void random(intAL& L,uint N,uint d,int lo,int hi);

template<class T,class S>
void getCarray(MT::List<T>& L,MT::Array<S>& Ca){
  Ca.resize(L.N);
  uint i;
  for(i=0;i<L.N;i++) Ca(i)=L.nodes(i)->p;
}

#ifdef MT_IMPLEMENTATION
#  include "list.cpp"
#endif

#endif
