/*  Copyright (C) 2000, 2006  Marc Toussaint (mtoussai@inf.ed.ac.uk)
    under the terms of the GNU LGPL (http://www.gnu.org/copyleft/lesser.html)
    see the `std.h' file for a full copyright statement  */

/*! \file array.h
    \brief core file: defines the MT::Array class */

#ifndef MT_array_h
#define MT_array_h

#include"std.h"
#undef max
#undef min

#define FOR1D(x,i) for(i=0;i<x.d0;i++)
#define FOR2D(x,i,j) for(i=0;i<x.d0;i++) for(j=0;j<x.d1;j++)

#define forAll(i,A) for(i=A.p;i!=A.pstop;i++)

extern int ARRAYOLDREAD;

//===========================================================================
//
// fwd declarations
//

namespace MT{
  template<class T> class Array;
}
template<class T> T product(const MT::Array<T>& v);
template<class T> void innerProduct(MT::Array<T>& x,const MT::Array<T>& y, const MT::Array<T>& z);
template<class T> MT::Array<T>& plusA(MT::Array<T>& x,const MT::Array<T>& y,const MT::Array<T>& z);

//===========================================================================
//
// Array class
//

namespace MT{
/*!\brief Simple array container to store arbitrary-dimensional arrays
    (tensors); can buffer more memory than necessary for faster
    resize; enables non-const reference of subarrays; enables fast
    memove for elementary types; implements many standard
    array/matrix/tensor operations. Please see the fully public attributes at the
    bottom of this page -- everthing is meant to be perfectly
    transparent. Interfacing with ordinary C-buffers is simple,
    e.g. via \c Array::referTo (Cbuffer,size) and \c Array::p and \c
    Array::pp. Please see also the reference for the \ref array.h
    header, which contains lots of functions that can be applied on
    Arrays. */

template<class T>
class Array{
public:
  typedef Array<T> ArrayT;

public:
  T *p;	    //!< the pointer on the linear memory allocated
  uint N;   //!< number of elements
  uint nd;  //!< number of dimensions
  uint d0;  //!< 0th dim
  uint d1;  //!< 1st dim
  uint d2;  //!< 2nd dim
  uint *d;  //!< pointer to dimensions (for nd<3 points to d0)
  T *pstop; //!< end of memory (pstop is already out of the bound)
  uint M;   //!< size of allocated memory (may be greater than N)
  bool reference;//!< true if this refers to some external memory
  T **pp;   //!< C-style 2D pointer (only valid if created with \c Array::getCarray!)
  Array<uint> *sparse; //!< to allocate element, column and row index lists (1+d0+d1 lists)
  bool temp;

  enum matrixType{ full, diag };
  matrixType mtype;
  static int sizeT;


public://!@name options

  /*!\brief be careful!!!; this option will modify the \c
      Array::resizeCopy : instead of calling the copy operator= for
      each element, the memory is moved with the \c memmove command -
      works only if you know what you do ... [default is true for
      simple types (like double, int, etc) and false for unknown
      ones] */
  bool memMove;

  /*!\brief if flexiMem is true (which is default!) the resize method will
      (1) at the first call allocate the exact amount of memory, (2)
      at further calls of increasing memory allocate twice the memory
      needed or (3) at further calls of decreasing memory only free
      the memory if the new size is smaller than a fourth */
  bool flexiMem;
  //@}


private:
  static char memMoveInit;
  void init(){
    reference=false;
    memMove=false;
    if(sizeT==-1) sizeT=sizeof(T);
    if(memMoveInit==-1){
      memMoveInit=0;
      if(typeid(T)==typeid(bool) ||
	typeid(T)==typeid(char) ||
	typeid(T)==typeid(unsigned char) ||
	typeid(T)==typeid(int) ||
	typeid(T)==typeid(unsigned int) ||
	typeid(T)==typeid(short) ||
	typeid(T)==typeid(unsigned short) ||
	typeid(T)==typeid(long) ||
	typeid(T)==typeid(unsigned long) ||
	typeid(T)==typeid(float) ||
	typeid(T)==typeid(double)) memMoveInit=1;
    }
    memMove=(memMoveInit==1);
    flexiMem=true;
    p=pstop=0;
    M=N=nd=d0=d1=d2=0;
    pp=0;
    sparse=0;
    d=&d0;
    temp=false;
    mtype=full;
  }


public://!@name constructors

  //! standard constructor -- this becomes an empty array
  Array(){ init(); }

  //! copy constructor
  Array(const ArrayT& a){ init(); operator=(a); }

  //! constructor with resize
  Array(uint i){ init(); resize(i); }

  //! constructor with resize
  Array(uint i,uint j){ init(); resize(i,j); }

  //! constructor with resize
  Array(uint i,uint j,uint k){ init(); resize(i,j,k); }

  //! constructor with resize
  Array(uint i,uint j,uint k,uint l){ init(); resize(4,&i); }

  //! this becomes a reference on a subdimension of \c a
  Array(const ArrayT& a,uint dim){ init(); referToSubDim(a,dim); }

  //! this becomes a reference on a subdimension of \c a
  Array(const ArrayT& a,uint i,uint j){ init(); referToSubDim(a,i,j); }

  //! this becomes a reference on the C-array \c p
  Array(T* p,uint size){ init(); referTo(p,size); }

  ~Array(){
    freeMEM();
  }
  

public://!@name resize

  //! frees all memory; this becomes an empty array
  void clear(){ freeMEM(); }

  //! resize 1D array, discard the previous contents
  ArrayT& resize(uint D0){ nd=1; d0=D0; resetD(); resizeMEM(d0,false); return *this; }

  //! resize but copy the previous contents
  ArrayT& resizeCopy(uint D0){ nd=1; d0=D0; resetD(); resizeMEM(d0,true); return *this; }

  //! resize and initialize with zero
  ArrayT& resizeZero(uint D0){ nd=1; d0=D0; resetD(); resizeMEM(d0,false); setZero(); return *this; }

  //! reshape the dimensionality (e.g. from 2D to 1D); throw an error if this actually requires to resize the memory
  ArrayT& reshape(uint D0){ CHECK(N==D0,"reshape must preserve total memory size"); nd=1; d0=D0; d1=d2=0; resetD(); return *this; }

  //! same for 2D ...
  ArrayT& resize(uint D0,uint D1){ nd=2; d0=D0; d1=D1; resetD(); resizeMEM(d0*d1,false); return *this; }

  //! ...
  ArrayT& resizeCopy(uint D0,uint D1){ nd=2; d0=D0; d1=D1; resetD(); resizeMEM(d0*d1,true); return *this; }

  //! ...
  ArrayT& reshape(uint D0,uint D1){ CHECK(N==D0*D1,"reshape must preserve total memory size"); nd=2; d0=D0; d1=D1; d2=0; resetD(); return *this; }

  //! same for 3D ...
  ArrayT& resize(uint D0,uint D1,uint D2){ nd=3; d0=D0; d1=D1; d2=D2; resetD(); resizeMEM(d0*d1*d2,false); return *this; }

  //! ...
  ArrayT& resizeCopy(uint D0,uint D1,uint D2){ nd=3; d0=D0; d1=D1; d2=D2; resetD(); resizeMEM(d0*d1*d2,true); return *this; }

  //! ...
  ArrayT& reshape(uint D0,uint D1,uint D2){ CHECK(N==D0*D1*D2,"reshape must preserve total memory size"); nd=3; d0=D0; d1=D1; d2=D2; resetD(); return *this; }


  //! resize to multi-dimensional tensor
  void resize(uint ND,uint *dim){
    nd=ND; d0=d1=d2=0; resetD();
    uint j,S;
    for(j=0;j<nd && j<3;j++){ (&d0)[j]=dim[j]; }
    if(nd>3){ d=new uint[nd];  memmove(d,dim,nd*sizeof(uint)); }
    for(S=1,j=0;j<nd;j++) S*=dim[j];
    resizeMEM(S,false);
  }

  void resize(const Array<uint> &newD){ resize(newD.N,newD.p); }

#ifndef MT_MSVC
  //! make it the same size as \c a
  void resizeAs(const ArrayT& a){
    nd=a.nd; d0=a.d0; d1=a.d1; d2=a.d2; resetD();
    if(nd>3){ d=new uint[nd];  memmove(d,a.d,nd*sizeof(uint)); }
    resizeMEM(a.N,false);
  }

  //! make it the same size as \c a and copy previous content
  void resizeCopyAs(const ArrayT& a){
    nd=a.nd; d0=a.d0; d1=a.d1; d2=a.d2; resetD();
    if(nd>3){ d=new uint[nd];  memmove(d,a.d,nd*sizeof(uint)); }
    resizeMEM(a.N,true);
  }
#endif

  //! resize to the same size as \c a
  template<class S>
  void resizeAs(const Array<S>& a){
    nd=a.nd; d0=a.d0; d1=a.d1; d2=a.d2; resetD();
    if(nd>3){ d=new uint[nd];  memmove(d,a.d,nd*sizeof(uint)); }
    resizeMEM(a.N,false);
  }

  //! make it the same size as \c a and copy previous content
  template<class S>
  void resizeCopyAs(const Array<S>& a){
    nd=a.nd; d0=a.d0; d1=a.d1; d2=a.d2; resetD();
    if(nd>3){ d=new uint[nd];  memmove(d,a.d,nd*sizeof(uint)); }
    resizeMEM(a.N,true);
  }

  //! return the size of memory allocated in bytes
  uint memsize(){ return M*sizeof(T); }


public://!@name tensor

  //! resize to multi-dimensional tensor
  void tensorResize(const Array<uint> &newD){ resize(newD.N,newD.p); }

  /*! makes X to be a distribution over the left leftmost-indexed
      variables and normalizes it */
  void tensorCondNormalize(int left){
    uint i,j,dl=1,dr;
    double sum;
    if(left>=0){ //normalize over the left variables
      for(j=0;j<(uint)left;j++) dl*=D(j);
      dr=N/dl;
      CHECK(dl*dr==N,"");
      for(i=0;i<dr;i++){
	sum=0.;
	for(j=0;j<dl;j++)  sum += p[j*dr + i];
	if(sum) for(j=0;j<dl;j++) p[j*dr + i] /= sum;
	else    for(j=0;j<dl;j++) p[j*dr + i] = 1./dl;
      }
    }else{ //normalize over the right variables
      for(j=0;j<(uint)-left;j++) dl*=D(j);
      dr=N/dl;
      CHECK(dl*dr==N,"");
      for(i=0;i<dl;i++){
	sum=0.;
	for(j=0;j<dr;j++)  sum += p[i*dr + j];
	if(sum) for(j=0;j<dr;j++) p[i*dr + j] /= sum;
	else    for(j=0;j<dr;j++) p[i*dr + j] = 1./dr;
      }
    }
  }

  /*! makes X to be a distribution over the left leftmost-indexed
      variables and normalizes it */
  void tensorCondMax(uint left){
    uint i,j,dl=1,dr,jmax;
    double pmax;
    for(j=0;j<left;j++) dl*=D(j);
    dr=N/dl;
    CHECK(dl*dr==N,"");
    for(i=0;i<dr;i++){
      jmax=0;
      pmax=p[jmax*dr + i];
      p[i]=0.;
      for(j=1;j<dl;j++){
	if(p[j*dr + i]>pmax){  jmax=j;  pmax=p[jmax*dr + i];  }
	p[j*dr + i]=0.;
      }
      p[jmax*dr + i]=1.;
    }
  }

  /*! makes X to be a distribution over the left leftmost-indexed
      variables and normalizes it */
  void tensorCondSoftMax(uint left,double beta){
    uint i;
    for(i=0;i<N;i++) p[i] = ::exp(beta*p[i]);
    tensorCondNormalize(left);
  }

  //! multi-dimensional (tensor) reference
  T& operator()(const Array<uint> &I) const{
    uint i,j;
    i=0;
    for(j=0;j<nd;j++) i = i*D(j) + I(j);
    return p[i];
  }

  //! I becomes the index tuple for the absolute index i
  void getIndexTuple(Array<uint> &I,uint i) const{
    uint j;
    CHECK(i<N,"out of range");
    I.resize(nd);
    I.setZero();
    for(j=nd;j--;){
      I.p[j] = i % d[j];
      i -= I.p[j];
      i /= d[j];
    }
  }

  //! return the k-th dimensionality
  uint D(uint k) const{
    CHECK(k<nd,"dimensionality range check error: " <<k <<"!<" <<nd);
    if(!d && k<3) return (&d0)[k]; else return d[k];
  }


public://!@name sparse arrays

  //! return fraction of non-zeros in the array
  double sparsity(){
    uint i,m=0;
    for(i=0;i<N;i++) if(elem(i)) m++;
    return ((double)m)/N;
  }

  //! make sparse: create the \ref sparse index
  void makeSparse(){
    CHECK(!sparse,"only once yet");
    uint n=0;
    if(nd==1){
      uint i;
      sparse=new Array<uint> [2];
      sparse[1].resize(d0); sparse[1]=-1;
      for(i=0;i<d0;i++) if(p[i]){
	sparse[0].append(i); //list of entries (maps n->i)
	sparse[1](i)=n;      //index list to entries (maps i->n)
	permute(i,n);
	n++;
      }
      N=n; resizeMEM(n,true);
      return;
    }
    if(nd==2){
      uint i,j;
      Array<uint> pair(2);
      sparse=new Array<uint> [1+d1+d0];
      for(i=0;i<d0;i++) for(j=0;j<d1;j++) if(p[i*d1+j]){
	pair(0)=i; pair(1)=j; sparse[0].append(pair);   sparse[0].reshape(n+1,2);
	permute(i*d1+j,n);
	//register entry in columns an row indices
	pair(0)=i; pair(1)=n; sparse[1+j]   .append(pair); sparse[1+j]   .reshape(sparse[1+j]   .N/2,2);
	pair(0)=j; pair(1)=n; sparse[1+d1+i].append(pair); sparse[1+d1+j].reshape(sparse[1+d1+j].N/2,2);
	n++;
      }
      N=n; resizeMEM(n,true);
      return;
    }
  }


public://!@name internal memory routines (probably not for external use)

  //! allocate memory (maybe using \ref flexiMem)
  void resizeMEM(uint n,bool copy){
    if(n==N) return;
    CHECK(!reference,"real resize of subarray is not allowed! (only a resize without changing memory size)");
    uint i;
    T *pold=p;
    if(true || !flexiMem){
      if(M!=n){
	if(n){
	  p=new T [n];
	  //p=(T*)malloc(n*sizeT);
	  if(!p) HALT("MT::Array failed memory allocation of "<<n*sizeT<<"bytes");
	  if(copy && !memMove) for(i=N<n?N:n;i--;) p[i]=pold[i];
	  if(copy && memMove) memmove(p,pold,sizeT*(N<n?N:n));
	}else p=0;
	CHECK((pold && M) || (!pold && !M),"");
	if(M) delete[] pold;
	//if(M) free(pold);
	N=M=n;
	pstop=p+N;
      }
    }else{
      if(n>0 && M==0){ //first time
	p=new T [n];
	if(!p) HALT("MT::Array failed memory allocation of "<<n*sizeT<<"bytes");
	M=n;
      }else if(n>M || 10+2*n<M/2){
	uint oversize=10+2*n;
	p=new T [oversize];
	if(!p) HALT("MT::Array failed memory allocation of "<<n*sizeT<<"bytes");
	if(copy && !memMove) for(i=N<n?N:n;i--;) p[i]=pold[i];
	if(copy && memMove) memmove(p,pold,sizeT*(N<n?N:n));
	if(M) delete[] pold;
	M=oversize;
      }
      N=n;
      pstop=p+N;
    }
    if(pp){ delete[] pp; pp=0; }
  }

  //! free all memory and reset all pointers and sizes
  void freeMEM(){
    if(M) delete[] p;
    //if(M) free(p);
    if(pp) delete[] pp;
    if(sparse) delete[] sparse;
    if(d && d!=&d0) delete[] d;
    p=pstop=0;
    M=N=nd=d0=d1=d2=0;
    sparse=0;
  }

  //! reset the dimensionality pointer d to point to &d0
  void resetD(){
    if(d && d!=&d0) delete[] d;
    d=&d0;
  }



public://!@name append, insert & remove

  //! append an (uninitialized) element to the array and return its reference -- the array becomes 1D!
  T& append(){
    if(nd==2 && d1==1)
      resizeCopy(d0+1,d1);
    else
      resizeCopy(N+1);
    return p[N-1];
  }

  //! append an element to the array -- the array becomes 1D!
  void append(const T& x){ append()=x; }

  //! append another array to the array (by copying it) -- the array might become 1D!
  void append(const ArrayT& x){
    uint oldN=N,i;
    if(nd==2 && d1==x.N) 
      resizeCopy(d0+1,d1);
    else
      resizeCopy(N+x.N);
    if(!memMove){
      for(i=0;i<x.N;i++) p[oldN+i]=x.p[i];
    }else memmove(p+oldN,x.p,sizeT*x.N);
  }

  //! append a C array to the array (by copying it) -- the array might become 1D!
  void append(const T *p,uint n){
    uint oldN=N,i;
    if(nd==2 && d1==n) 
      resizeCopy(d0+1,d1);
    else
      resizeCopy(N+n);
    if(!memMove){
      for(i=0;i<n;i++) p[n+i]=p[i];
    }else memmove(p+oldN,p,sizeT*n);
  }

  //! the array contains `copies' copies of the old one
  void replicate(uint copies){
    if(copies<2) return;
    uint i,oldN=N;
    resizeCopy(copies*N);
    if(!memMove){
      HALT("not implemented yet");
    }else{
      for(i=0;i<copies;i++) memmove(p+i*oldN,p,sizeT*oldN);
    }
  }

  //! inserts x at the position i -- the array becomes 1D! [only with memMove!]
  void insert(uint i,const T& x){
    CHECK(memMove,"only with memMove");
    uint Nold=N;
    resizeCopy(Nold+1);
    if(i<Nold) memmove(p+i+1,p+i,sizeT*(Nold-i));
    p[i]=x;
  }

  //! remove (delete) a subsequence of the array -- the array becomes 1D!  [only with memMove!]
  void remove(uint i,uint n=1){
    CHECK(memMove,"only with memMove");
    if(N>i+n) memmove(p+i,p+i+n,sizeT*(N-i-n));
    resizeCopy(N-n);
  }

  //! remove some element by permuting the last element in its place! -- the array becomes 1D!
  void removePerm(uint i){
    p[i]=p[N-1];
    resizeCopy(N-1);
  }

  //! replace n elements at pos i by the sequence x -- the array becomes 1D!  [only with memMove!]
  void replace(uint i,uint n,const ArrayT& x){
    CHECK(memMove,"only with memMove");
    uint Nold=N;
    if(n==x.N){
      memmove(p+i,x.p,sizeT*(x.N));
    }else if(n>x.N){
      memmove(p+i+x.N,p+i+n,sizeT*(Nold-i-n));
      if(i+n<Nold) memmove(p+i,x.p,sizeT*(x.N));
      resizeCopy(Nold-n+x.N);
    }else{
      resizeCopy(Nold+x.N-n);
      if(i+n<Nold) memmove(p+i+x.N,p+i+n,sizeT*(Nold-i-n));
      memmove(p+i,x.p,sizeT*(x.N));
    }
  }

  //! deletes the i-th row [must be 2D]
  void delRow(uint i){
    CHECK(memMove,"only with memMove");
    CHECK(nd==2,"only for matricies");
    CHECK(i<d0,"range check error");
    uint n=d1;
    if(i+1<d0) memmove(p+i*n,p+(i+1)*n,sizeT*(d0-i-1)*n);
    resizeCopy(d0-1,n);
  }

  //! deletes the i-th column [must be 2D]
  void delColumns(uint i,uint k=1){
    CHECK(memMove,"only with memMove");
    CHECK(nd==2,"only for matricies");
    CHECK(i<d1,"range check error");
    uint n=d1;
    for(uint j=0;j<d0;j++){
      memmove(p+j*(n-k)  ,p+j*n      ,sizeT*i);
      memmove(p+j*(n-k)+i,p+j*n+(i+k),sizeT*(n-i-k));
    }
    resizeCopy(d0,n-k);
  }

  //! inserts k columns at the i-th column [must be 2D]
  void insColumns(uint i,uint k=1){
    CHECK(memMove,"only with memMove");
    CHECK(nd==2,"only for matricies");
    CHECK(i<=d1,"range check error");
    uint n=d1;
    resizeCopy(d0,n+k);
    for(uint j=d0;j--;){
      memmove(p+j*d1+(i+k),p+j*n+i,sizeT*(n-i));
      memset (p+j*d1+i    ,0      ,sizeT*k);
      memmove(p+j*d1      ,p+j*n  ,sizeT*i);
    }
  }
  
  //! inserts k rows at the i-th row [must be 2D]
  void insRows(uint i,uint k=1){
    CHECK(memMove,"only with memMove");
    CHECK(nd==2,"only for matricies");
    CHECK(i<=d0,"range check error");
    uint n=d0;
    resizeCopy(d0+k,d1);
    memmove(p+(i+k)*d1,p+i*d1,sizeT*d1*(n-i));
    memset (p+ i   *d1,0     ,sizeT*d1*k);
  }

public://!@name access operations

  //! the \c ith element
  T& elem(uint i) const{
    CHECK(i<N,"range error ("<<i<<">="<<N<<")");
    return p[i];
  }

  /* scalar reference (legal iff N==1) */
  /*operator T&() const{ 
    CHECK(N==1,"scalar reference ("<<N<<"!=1)");
    return *p;
  }*/

  //! 1D reference
  T& operator()(uint i) const{ 
    CHECK((nd==1 || nd>3) && i<N,
	  "1D range error ("<<nd<<"=1,"<<i<<"<"<<d0<<")");
    return p[i];
  }

  //! 2D reference
  T& operator()(uint i,uint j) const{
    CHECK(nd==2 && i<d0 && j<d1 && !sparse,
	  "2D range error ("<<nd<<"=2,"<<i<<"<"<<d0<<","<<j<<"<"<<d1<<")");
    return p[i*d1+j];
  }

  //! 3D reference
  T& operator()(uint i,uint j,uint k) const{
    CHECK(nd==3 && i<d0 && j<d1 && k<d2 && !sparse,
	  "3D range error ("<<nd<<"=3,"<<i<<"<"<<d0<<","<<j<<"<"<<d1<<","<<k<<"<"<<d2<<")");
    return p[(i*d1+j)*d2+k];
  }

  //! get a subarray; use in conjuction with operator()() to get a reference
  ArrayT operator[](uint i) const{ return Array(*this,i); }

  //! convert a subarray into a reference (e.g. a[3]()+=.123)
  ArrayT& operator()(){ return (*this); }

  //! reference to the max entry
  T& max() const{ uint i,m=0; for(i=0;i<N;i++) if(p[i]>p[m]) m=i; return p[m]; }

  //! reference to the min entry
  T& min() const{ uint i,m=0; for(i=0;i<N;i++) if(p[i]<p[m]) m=i; return p[m]; }

  /*!\brief the index of the maxium; precondition: the comparision operator
      > exists for type T */ 
  uint maxIndex(){ uint i,m=0; for(i=0;i<N;i++) if(p[i]>p[m]) m=i; return m; }

  /*!\brief the index of the maxium; precondition: the comparision operator
      > exists for type T */ 
  void maxIndex(uint& i,uint& j){ CHECK(nd==2,"needs 2D array"); j=maxIndex(); i=j/d1; j=j%d1; }

  /*!\brief the index of the maxium; precondition: the comparision operator
      > exists for type T */ 
  void maxIndex(uint& i,uint& j,uint& k){ CHECK(nd==3,"needs 3D array"); k=maxIndex(); i=k/(d1*d2); k=k%(d1*d2); j=k/d2; k=k%d2; }

  //! ...
  uint minIndex(){ uint i,m=0; for(i=0;i<N;i++) if(p[i]<p[m]) m=i; return m; }

  int contains(T& x){ uint i; for(i=0;i<N;i++) if(p[i]==x) return i; return -1; }

  /*!\brief a sub array of a 1D Array (corresponds to matlab [i:I]); when
      the upper limit I is -1, it is replaced by the max limit (like
      [i:]) */
  ArrayT sub(uint i,int I) const{
    CHECK(nd==1,"1D range error ");
    ArrayT x;
    if(I==-1) I=d0-1;
    CHECK(i<=(uint)I,"lower limit higher than upper!");
    x.resize(I-i+1);
    uint k;
    for(k=i;k<=(uint)I;k++) x(k-i)=operator()(k);
    return x;
  }

  /*!\brief copies a sub array of a 2D Array (corresponds to matlab [i:I,j:J]);
      when the upper limits I or J are -1, they are replaced by the
      max limit (like [i:,j:]) */
  ArrayT sub(uint i,int I,uint j,int J) const{
    CHECK(nd==2,"2D range error ");
    ArrayT x;
    if(I==-1) I=d0-1;
    if(J==-1) J=d1-1;
    CHECK(i<=(uint)I && j<=(uint)J,"lower limit must be higher than upper!");
    x.resize(I-i+1,J-j+1);
    uint k,l;
    for(k=i;k<=(uint)I;k++) for(l=j;l<=(uint)J;l++) x(k-i,l-j)=operator()(k,l);
    return x;
  }

  /*!\brief copies a sub array of a 2D Array (corresponds to matlab [i:I,j:J]);
      when the upper limits I or J are -1, they are replaced by the
      max limit (like [i:,j:]) */
  ArrayT sub(uint i,int I,uint j,int J,uint k,int K) const{
    CHECK(nd==3,"3D range error ");
    ArrayT x;
    if(I==-1) I=d0-1;
    if(J==-1) J=d1-1;
    if(K==-1) K=d2-1;
    CHECK(i<=(uint)I && j<=(uint)J && k<=(uint)K,"lower limit must be higher than upper!");
    x.resize(I-i+1,J-j+1,K-k+1);
    uint ii,jj,kk;
    for(ii=i;ii<=(uint)I;ii++) for(jj=j;jj<=(uint)J;jj++)  for(kk=k;kk<=(uint)K;kk++) x(ii-i,jj-j,kk-k)=operator()(ii,jj,kk);
    return x;
  }

  /*!\brief copies a selection of columns from a 2D array, the first index (rows)
    runs from i to I (as explained above) while the second index runs
    over the columns explicitly referred to by cols. (col doesn't have
    to be ordered or could also contain some columns multiply) */
  ArrayT sub(uint i,int I,Array<uint> cols) const{
    CHECK(nd==2,"2D range error ");
    ArrayT x;
    if(I==-1) I=d0-1;
    CHECK(i<=(uint)I,"lower limit higher than upper!");
    x.resize(I-i+1,cols.N);
    uint k,l;
    for(k=i;k<=(uint)I;k++) for(l=0;l<cols.N;l++) x(k-i,l)=operator()(k,cols(l));
    return x;
  }


public://!@name C-array interfacing

  //! allocates, sets and returns the \c Array::pp pointer (of type \c T**)
  T** getCarray() const{
    CHECK(nd>=2,"only 2D or higher-D arrays gives C-array of type T**");
    if(pp) return pp;
    ((ArrayT*)this)->pp=new T* [d0];
    uint skip;
    if(nd==2) skip=d1; else skip=N/d0;
    for(uint i=0;i<d0;i++) pp[i]=p+i*skip;
    return pp;
  }

  //! makes this array a reference to the C buffer
  void referTo(const T *buffer,uint n){
    freeMEM();
    reference=true;
    nd=1; d0=n; d1=d2=0; N=n;
    p=(T*)buffer;
    pstop=p+N;
  }

#if 0
  /*!\brief returns an ordinary 2-dimensional C-pointer to the Array content.
      Requires the Array<T*> as buffer. */
  T** getPointers(Array<T*>& array2d) const{
    CHECK(nd==2,"only 2D array gives C-array of type T**");
    array2d.resize(d0);
    for(uint i=0;i<d0;i++) array2d(i)=p+i*d1;
    return array2d.p;
  }

  //! returns an ordinary 3-dimensional C-pointer-array
  T*** getPointers(Array<T**>& array3d,Array<T*>& array2d) const{
    CHECK(nd==3,"only 3D array gives C-array of type T*** ");
    array2d.resize(d0,d1);
    for(uint i=0;i<d0;i++){
      for(uint j=0;j<d1;j++) array2d(i,j)=&operator()(i,j,0);
      array3d(i)=&array2d(i,0);
    }
    return array3d.p;
  }
#endif

  //! copies the ordinary C buffer into this array
  template<class S>
  void copy(S *buffer,uint D0){
    resize(D0);
    uint i;
    if(memMove && typeid(T)==typeid(S))
      memmove(p,buffer,sizeT*d0);
    else for(i=0;i<d0;i++) operator()(i)=(T)buffer[i];
  }

  //! copies the ordinary C buffer into this array
  template<class S>
  void copy(S **buffer,uint D0,uint D1){
    resize(D0,D1);
    uint i,j;
    for(i=0;i<d0;i++){
      if(memMove && typeid(T)==typeid(S))
	memmove(p+i*d1,buffer[i],sizeT*d1);
      else for(j=0;j<d1;j++) operator()(i,j)=(T)buffer[i][j];
    }
  }

  /*!\brief copies this array into an ordinary C buffer !!be sure that
      buffer has the right size!! */
  template<class S>
  void copyInto(S *buffer){
    CHECK(nd==1,"can only copy 1D Array into 1D C-array");
    uint i;
    if(memMove && typeid(T)==typeid(S)) memmove(buffer,p,sizeT*d0);
    else for(i=0;i<d0;i++) buffer[i]=(S)operator()(i);
  }

  /*!\brief copies this array into an ordinary C buffer !!be sure that
      buffer has the right size!! */
  template<class S>
  void copyInto2D(S **buffer){
    CHECK(nd==2,"can only copy 2D Array into 2D C-array");
    uint i,j;
    for(i=0;i<d0;i++){
      if(memMove && typeid(T)==typeid(S)) memmove(buffer[i],p+i*d1,sizeT*d1);
      else for(j=0;j<d1;j++) buffer[i][j]=(S)operator()(i,j);
    }
  }


public://!@name assignments

  //! set all elements to value \c v 
  ArrayT& operator=(const T v){
    uint i;
    //if(memMove && typeid(T)==typeid(T)) memset(p,*((int*)&v),N); else
    for(i=0;i<N;i++) p[i]=v;
    return *this;
  }

  //! copy operator
  ArrayT& operator=(const ArrayT& a){
    //if(a.temp){ takeOver(*((ArrayT*)&a)); return *this; }
    resizeAs(a);
    uint i;
    if(memMove) memmove(p,a.p,sizeT*N);
    else for(i=0;i<N;i++) p[i]=a.p[i];
    return *this;
  }

  //! reads it from the C string (by using std::istringstream)
  void set(const char* str){
    std::istringstream s(str);
    readRaw(s);
  }

  //! generic copy; precondition: operator (T)=(S) exists
  template<class S>
  ArrayT& copy(const Array<S>& a){
    resizeAs(a);
    uint i;
    if(memMove && typeid(T)==typeid(S)) memmove(p,a.p,sizeT*N);
    else for(i=0;i<N;i++) p[i]=(T)a.p[i];
    return *this;
  }

  /*!\brief same as memset(p,0,sizeT*N); precondition: memMove is
      true! */
  void setZero(byte zero=0){
    CHECK(memMove,"can set array's memory to zero only if memMove option is true");
    memset(p,zero,sizeT*N);
  }

  //! set all entries to same value x [default: don't change dimension]
  void setUni(const T& x,int d=-1){
    if(d!=-1) resize(d);
    uint i;
    for(i=0;i<N;i++) operator()(i)=x;
  }

  /*!\brief becomes the n-dim identity matrix [default:
      don't change dimension (already has to be squared matrix)] */
  void setId(int d=-1){
    CHECK(d!=-1 || (nd==2 && d0==d1),"need squared matrix to set to identity");
    if(d!=-1) resize(d,d);
    setZero();
    for(uint i=0;i<d0;i++) operator()(i,i)=1.;
  }

  void setDiag(const T& x,int d=-1){
    CHECK(d!=-1 || (nd==2 && d0==d1),"need squared matrix to set to identity");
    if(d!=-1) resize(d,d);
    setZero();
    uint i;
    for(i=0;i<d0;i++) operator()(i,i)=x;
  }

  void setBlock(const ArrayT& A,const ArrayT& B,const ArrayT& C,const ArrayT& D){
    CHECK(A.nd==2 && B.nd==2 && C.nd==2 && D.nd==2,"");
    CHECK(A.d0==B.d0 && A.d1==C.d1 && B.d1==D.d1 && C.d0==D.d0,"");
    uint i,j,a=A.d0,b=A.d1;
    resize(A.d0+C.d0,A.d1+B.d1);
    for(i=0;i<A.d0;i++) for(j=0;j<A.d1;j++) operator()(i  ,j  )=A(i,j);
    for(i=0;i<B.d0;i++) for(j=0;j<B.d1;j++) operator()(i  ,j+b)=B(i,j);
    for(i=0;i<C.d0;i++) for(j=0;j<C.d1;j++) operator()(i+a,j  )=C(i,j);
    for(i=0;i<D.d0;i++) for(j=0;j<D.d1;j++) operator()(i+a,j+b)=D(i,j);
  }

  void setBlock(const ArrayT& B,uint lo0,uint lo1){
    CHECK(nd==2 && B.nd==2 && lo0+B.d0<=d0 && lo1+B.d1<=d1,"");
    uint i,j;
    for(i=0;i<B.d0;i++) for(j=0;j<B.d1;j++) operator()(lo0+i,lo1+j)=B(i,j);
  }

  void getBlock(ArrayT& B,uint lo0,uint lo1) const{
    CHECK(nd==2 && B.nd==2 && lo0+B.d0<=d0 && lo1+B.d1<=d1,"");
    uint i,j;
    for(i=0;i<B.d0;i++) for(j=0;j<B.d1;j++) B(i,j)=operator()(lo0+i,lo1+j);
  }

  void setBlock(const ArrayT& B,uint lo){
    CHECK(nd==1 && B.nd==1 && lo+B.N<=N,"");
    uint i;
    for(i=0;i<B.N;i++) operator()(lo+i)=B(i);
  }

  void getBlock(ArrayT& B,uint lo) const{
    CHECK(nd==1 && B.nd==1 && lo+B.N<=N,"");
    uint i;
    for(i=0;i<B.N;i++) B(i)=operator()(lo+i);
  }

  void setStraightPerm(int n=-1){
    if(n!=-1) resize(n);
    for(uint i=0;i<N;i++) elem(i)=i;
  }

  //! make this array a reference to the array \c a
  void referTo(const ArrayT& a){
    freeMEM();
    reference=true; memMove=a.memMove; flexiMem=a.flexiMem;
    N=a.N; nd=a.nd; d0=a.d0; d1=a.d1; d2=a.d2;
    p=a.p;
    pstop=a.pstop;
  }

  //! make this array a subarray reference to \c a
  void referToSubRange(const ArrayT& a,uint i,int I){
    CHECK(a.nd==1,"not implemented yet");
    freeMEM();
    reference=true; memMove=a.memMove; flexiMem=a.flexiMem;
    if(I==-1) I=a.N-1;
    N=I+1-i; nd=1; d0=N; d1=0; d2=0;
    p=a.p+i;
    pstop=p+N;
  }

  //! make this array a subarray reference to \c a
  void referToSubDim(const ArrayT& a,uint dim){
    CHECK(a.nd>1,"can't create subarray of array less than 2 dimensions");
    freeMEM();
    reference=true; memMove=a.memMove; flexiMem=a.flexiMem;
    if(a.nd==2){
      nd=1; d0=a.d1; d1=d2=0; N=d0;
      p=&a(dim,0);
      pstop=p+N;
    }
    if(a.nd==3){
      nd=2; d0=a.d1; d1=a.d2; d2=0; N=d0*d1;
      p=&a(dim,0,0);
      pstop=p+N;
    }
  }

  //! make this array a subarray reference to \c a
  void referToSubDim(const ArrayT& a,uint i,uint j){
    CHECK(a.nd>2,"can't create subsubarray of array less than 3 dimensions");
    freeMEM();
    reference=true; memMove=a.memMove; flexiMem=a.flexiMem;
    if(a.nd==3){
      nd=1; d0=a.d2; d1=0; d2=0; N=d0;
      p=&a(i,j,0);
      pstop=p+N;
    }
  }

  //! redirects a reference (previously constructed) to a new sub dim of a
  void redirect(const ArrayT& a,uint i){
    CHECK(reference && a.nd==2,"can only do that hack with references");
    p=a.p+i*a.d1;
  }

  /*! takes over the memory buffer from a; afterwards, this is a
    proper array with own memory and a is only a reference on the
    memory */
  void takeOver(ArrayT& a){
    freeMEM();
    memMove=a.memMove; flexiMem=a.flexiMem;
    N=a.N; nd=a.nd; d0=a.d0; d1=a.d1; d2=a.d2;
    p=a.p;
    pstop=a.pstop;
    M=a.M;
    a.reference=true;
    a.M=0;
  }


public://!@name stuff for generic looping
  void First(T*& i) const{ i=p; }
  void Last(T*& i) const{ i=(p+N)-1; }
  void Next(T*& i) const{ i++; }
  bool CheckNext(T*& i) const{ return i!=pstop; }
  void Prev(T*& i) const{ i--; }
  bool CheckPrev(T*& i) const{ return i!=p-1; }


public://!@name permutations

  //! permute the elements \c i and \c j
  void permute(uint i,uint j){ T x=p[i]; p[i]=p[j]; p[j]=x; }

  //! permute all elements randomly
  void randomPermute(){
    int j,r;
    for(j=N-1;j>=1;j--){
      r=rnd(j+1);
      permute(r,j);
    }
  }
  
  //! push all elements backward
  void push(int offset,bool circle=true){
    static ArrayT tmp;
    CHECK(memMove,"pushing only works with memMove enabled");
    uint m=abs(offset);
    if(circle){
      tmp.resize(m);
      if(offset>0){
        memmove(tmp.p,p+N-m,sizeT*m);
        memmove(p+m,p,sizeT*(N-m));
        memmove(p,tmp.p,sizeT*m);
      }
      if(offset<0){
        memmove(tmp.p,p,sizeT*m);
        memmove(p,p+m,sizeT*(N-m));
        memmove(p+(N-m),tmp.p,sizeT*m);
      }
    }
  }


public://!@name I/O

  //! prototype for operator<<
  void write(std::ostream& os,char *ELEMSEP=" ",char *LINESEP="\n ",bool binary=false) const{
    CHECK(!binary || memMove,"binary write works only for memMoveable data");
    uint i,j,k;
    if(!IOraw){//write tag
      if(binary){
	if(nd==0){ os <<"<0>"; return; }
	if(nd==1)  os <<"<1:" <<d0;
	if(nd==2)  os <<"<2:" <<d0 <<',' <<d1;
	if(nd==3)  os <<"<3:" <<d0 <<',' <<d1 <<',' <<d2;
      }else{
	if(nd==0){ os <<"[]"; return; }
	os <<"[";
      }
    }
    //write data
    if(binary){
      os.write((char*)p,sizeT*N);
    }else{
      //os <<std::setprecision(4) <<std::fixed <<std::setfill(' ');
      if(nd==1){
	//os <<' ';
	for(i=0;i<N;i++) os <<ELEMSEP <<std::setw(4) <<operator()(i);
      }
      if(nd==2) for(j=0;j<d0;j++){
	if(j) os <<LINESEP;
	for(i=0;i<d1;i++) os <<ELEMSEP <<std::setw(4) <<operator()(j,i);
      }
      if(nd==3) for(k=0;k<d0;k++){
	if(k) os <<LINESEP;
	for(j=0;j<d1;j++){
	  for(i=0;i<d2;i++) os <<ELEMSEP <<std::setw(4) <<operator()(k,j,i);
	  os <<LINESEP;
	}
      }
      if(nd>3){
	CHECK(d && d!=&d0,"")
	Array<uint> I;
	for(i=0;i<N;i++){
	  if(i && !(i%d[nd-1])) os <<LINESEP;
	  if(nd>1 && !(i%(d[nd-2]*d[nd-1]))){
	    getIndexTuple(I,i);
	    os <<LINESEP <<'<' <<i <<':' <<I(0);
	    for(j=1;j<nd;j++) os <<' ' <<I(j);
	    os <<'>' <<LINESEP;
	  }
	  os <<ELEMSEP <<std::setw(7) <<operator()(i);
	}
      }
    }
    if(!IOraw){//write end tag
      if(binary){
	os <<">";
      }else{
	os <<" ]";
      }
      //if(nd>1) os <<std::endl;
    }
  }

  //! prototype for operator>>
  void read(std::istream& is,bool binary=false){
    CHECK(!binary || memMove,"binary write works only for memMoveable data");
    if(ARRAYOLDREAD){ readOld(is); return; }
    uint d,i,j,k;
    char c;
    T x;
    MT::skip(is);
    is.get(c);
    switch(c){
    case '<':
      is >>d;
      if(is.fail()) HALT ("could not read array tag");
      if(d==0){ is >>">"; return; }
      if(d==1){
	is >>":" >>i;
	if(is.fail()) HALT ("could not read array's dimensions");
	resize(i);
      }
      if(d==2){
	is >>":" >>i >>"," >>j;
	if(is.fail()) HALT ("could not read array's dimensions");
	resize(i,j);
      }
      if(d==3){
	is >>":" >>i >>"," >>j >>"," >>k;
	if(is.fail()) HALT ("could not read array's dimensions");
	resize(i,j,k);
      }
      //read data
      if(binary){
	is.read((char*)p,sizeT*N);
      }else{
	for(i=0;i<N;i++){
	  if(is.fail()) HALT("could not read "<<i<<"-th element of an array");
	  is >>p[i];
	}
      }
      is >>">";
      if(is.fail()) HALT ("could not read array end tag");
      break;
    case '[':
    default:
      if(c!='[') is.putback(c);
      uint i=0;
      d=0;
      for(;;){
	MT::skip(is," \r\t");
	is.get(c);
	if(c==']' || !is.good()){ is.clear(); break; }
	if(c==';' || c=='\n'){ //set an array width
	  if(!d) d=i; else CHECK(!(i%d),"Array::read: mis-structured array in row "<<i/d);
	  continue;
	}
	if(c!=',') is.putback(c);
	is >>x;
	if(!is.good()){ is.clear(); break; }
	if(i>=N) resizeCopy(i+1000);
	elem(i)=x;
	i++;
      }
      resizeCopy(i);
      if(d){
	CHECK(!(N%d),"Array::read: mis-structured array in last row"); 
	reshape(N/d,d);
      }
      break;
    }
  }

  void readOld(std::istream& is){
    uint d,i,j,k;
    char c;
    MT::skip(is);
    is.get(c);
    switch(c){
    case '[':
      is >>d;
      if(is.fail()) HALT ("could not read array tag");
      if(d==0){ is >>">"; return; }
      if(d==1){
	is >>":" >>i;
	if(is.fail()) HALT ("could not read array's dimensions");
	resize(i);
      }
      if(d==2){
	is >>":" >>i >>"," >>j;
	if(is.fail()) HALT ("could not read array's dimensions");
	resize(i,j);
      }
      if(d==3){
	is >>":" >>i >>"," >>j >>"," >>k;
	if(is.fail()) HALT ("could not read array's dimensions");
	resize(i,j,k);
      }
      is >>"]";
      if(is.fail()) HALT ("could not read array end tag");
      for(i=0;i<N;i++){
	if(is.fail()) HALT("could not read "<<i<<"-th element of an array");
	is >>p[i];
      }
      break;
    default:
      NIY;
    }
  }

  void writeDim(std::ostream& os) const{
    uint i;
    os <<'<';
    if(nd) os <<D(0); else os <<0;
    for(i=1;i<nd;i++) os <<' ' <<D(i);
    os <<'>';
  }

  void readDim(std::istream& is){
    char c;
    uint ND,dim[10];
    is >>"<";
    for(ND=0;;ND++){
      is >>dim[ND];
      is.get(c);
      if(c=='>') break;
      CHECK(c==' ',"error in reading dimensionality");
    }
    resize(ND+1,dim);
  }

  //! IO-manip: arrays will be streamed as raw (without tags)
  const ArrayT& ioraw() const{ IOraw=true; return *this; }

  //! IO-manip: arrays will be streamed non-raw (with tags)
  const ArrayT& ionoraw() const{ IOraw=false; return *this; }

  //! read array from stream without knowing yet the dimensionalities
  void readRaw(std::istream& is,uint d1=0,uint d2=0){
    resize(0);
    T x;
    for(;;){
      MT::skip(is);
      if(!is.good()) break;
      is >>x;
      if(is.fail()) HALT("error when reading "<<N<<"th element in readRaw");
      if(d1 && d2 && N>=d1*d2) break;
      append(x);
    }
    if(d1 && N>d1){
      CHECK(!(N%d1),"read "<<N<<"elements - not dividable by "<<d1);
      resizeCopy(N/d1,d1);
    }
  }
  void readRaw(const char* filename,uint d1=0,uint d2=0){
    std::ifstream is; MT::open(is,filename); readRaw(is,d1,d2);
  }
  void writeRaw(const char* filename){
    std::ofstream os; MT::open(os,filename); os <<(*this).ioraw();
  }
  //@}

  //! transpose
  friend ArrayT operator~(const ArrayT& y){ ArrayT x; transpose(x,y); x.temp=true; return x; }
  //! negative
  friend ArrayT operator-(const ArrayT& y){ ArrayT x; negative(x,y);  x.temp=true; return x; }
  //! outer product (notation borrowed from the wedge product, though not anti-symmetric)
  friend ArrayT operator^(const ArrayT& y,const ArrayT& z){ ArrayT x; outerProduct(x,y,z); x.temp=true; return x; }

  //! inner product
  friend ArrayT operator*(const ArrayT& y,const ArrayT& z){ ArrayT x; innerProduct(x,y,z); x.temp=true; return x; }
  //! scalar multiplication
  friend ArrayT operator*(const ArrayT& y,T z){             ArrayT x; multS(x,y,z); x.temp=true; return x; }
  //! scalar multiplication
  friend ArrayT operator*(T y,const ArrayT& z){             ArrayT x; multS(x,y,z); x.temp=true; return x; }


#define BinaryOperator( op, name) \
  friend ArrayT operator op(const ArrayT& y,const ArrayT& z){ ArrayT x; name(x,y,z); x.temp=true; return x; } \
  friend ArrayT operator op(T y,const ArrayT& z){             ArrayT x; name##S(x,y,z); x.temp=true; return x; } \
  friend ArrayT operator op(const ArrayT& y,T z){             ArrayT x; name##S(x,y,z); x.temp=true; return x; } \

  BinaryOperator( + , plusA)
  BinaryOperator( - , minusA)
  BinaryOperator( % , mult)
  BinaryOperator( / , div)
#undef BinaryOperator

#define CompoundAssignmentOperator( op )                                  \
  friend ArrayT& operator op (ArrayT& x,const MT::Array<T>& y){           \
    CHECK(x.N==y.N,							  \
    "binary operator on different array dimensions ("<<x.N<<", "<<y.N<<")");\
    T *xp=x.p;								  \
    const T *yp=y.p;							  \
    for(; xp!=x.pstop; xp++,yp++) *xp op *yp;                             \
    x.temp=true; return x;                                                \
  }                                                                       \
  friend ArrayT& operator op ( ArrayT& x, T y ){                          \
    T *xp=x.p;								  \
    for(; xp!=x.pstop; xp++) *xp op y;                                    \
    x.temp=true; return x;                                                \
  }                                                                       \

  CompoundAssignmentOperator( |= )
  CompoundAssignmentOperator( ^= )
  CompoundAssignmentOperator( &= )
  CompoundAssignmentOperator( += )
  CompoundAssignmentOperator( -= )
  CompoundAssignmentOperator( *= )
  CompoundAssignmentOperator( /= )
  CompoundAssignmentOperator( %= )
#undef CompoundAssignmentOperator
};
}


//===========================================================================
//
//!@name standard types
//

typedef MT::Array<double> arr;
typedef MT::Array<double> doubleA;
typedef MT::Array<float>  floatA;
typedef MT::Array<uint>   uintA;
typedef MT::Array<int>    intA;
typedef MT::Array<char>   charA;
typedef MT::Array<byte>   byteA;
typedef MT::Array<bool>   boolA;
typedef MT::Array<MT::String> stringA;


//===========================================================================
//
//!@name basic operators
//

//! calls MT::Array<T>::read
template<class T> inline std::istream& operator>>(std::istream& is,MT::Array<T>& x){ x.read(is);return is; }

//! allows a notation such as x <<"[0 1; 2 3]"; to initialize an array x 
template<class T> inline MT::Array<T>& operator<<(MT::Array<T>& x,const char* str){ std::istringstream ss(str); ss >>x; return x; }

//! calls MT::Array<T>::write
template<class T> inline std::ostream& operator<<(std::ostream& os,const MT::Array<T>& x){
  x.write(os); return os;
}

//! check whether this and \c a have same dimensions
template<class T,class S>
inline bool samedim(const MT::Array<T>& a,const MT::Array<S>& b){
  return (b.nd==a.nd && b.d0==a.d0 && b.d1==a.d1 && b.d2==a.d2);
}

//! check for Nans in the array (checks x.elem(i)==x.elem(i) for all elements)
template<class T>
void checkNan(const MT::Array<T>& x){
  for(uint i=0;i<x.N;i++){
    //CHECK(x.elem(i)!=NAN,"found a NaN" <<x.elem(i) <<'[' <<i <<']');
    CHECK(x.elem(i)==x.elem(i),"inconsistent number: " <<x.elem(i) <<'[' <<i <<']');
  }
}

//! equal in size and all elements
template<class T>
inline bool operator==(const MT::Array<T>& v,const MT::Array<T>& w){
  if(!samedim(v,w)) return false;
  const T *iv,*iw;
  for(iv=v.p,iw=w.p; iv!=v.pstop; iv++,iw++)
    if (*iv != *iw) return false;
  return true;
}

//! equal in size and all elements
template<class T>
inline bool operator==(const MT::Array<T>& v,const T *w){
  const T *iv,*iw;
  for(iv=v.p,iw=w; iv!=v.pstop; iv++,iw++)
    if (*iv != *iw) return false;
  return true;
}

//! not equal
template<class T>
inline bool operator!=(const MT::Array<T>& v,const MT::Array<T>& w){
  return !(v==w);
}

//! lexical comparison
template<class T>
inline bool operator<(const MT::Array<T>& v,const MT::Array<T>& w){
  if(v.N==w.N){
    for(uint i=0;i<v.N;i++){
      if(v.p[i]>w.p[i]) return false;
      if(v.p[i]<w.p[i]) return true;
    }
    return false; //they are equal
  }
  return v.N<w.N;
}


//===========================================================================
//
//!@name arithmetic operators
//


//---------- unary operators

#define UnaryOperation( name, op )                                        \
template<class T>                                                         \
inline MT::Array<T>& name (MT::Array<T>& x,const MT::Array<T>& y){        \
  x.resizeAs(y);							  \
  T *xp=x.p;                                                              \
  const T *iy=y.p;                                                        \
  for(; iy!=y.pstop; iy++,xp++) *xp= op *iy;			  	  \
  return x;                                                               \
}

UnaryOperation( negative, - )
#undef UnaryOperator

//---------- binary function

#define BinaryOperation( name, op )                                       \
template<class T>                                                         \
inline MT::Array<T>& name(MT::Array<T>& x,const MT::Array<T>& y,const MT::Array<T>& z){ \
  CHECK(y.N==z.N,							  \
  "binary operator on different array dimensions ("<<y.N<<", "<<z.N<<")");\
  x.resizeAs(y);							  \
  T *xp=x.p;								  \
  const T *zp=z.p,*yp=y.p;						  \
  for(; xp!=x.pstop; xp++,yp++,zp++) *xp = *yp op *zp;                    \
  return x;                                                               \
}                                                                         \
                                                                          \
template<class T>                                                         \
inline MT::Array<T>& name##S(MT::Array<T>& x,const MT::Array<T>& y,T z){     \
  x.resizeAs(y);							  \
  T *xp=x.p;								  \
  const T *yp=y.p;							  \
  for(; xp!=x.pstop; xp++,yp++) *xp = *yp op z;                           \
  return x;                                                               \
}                                                                         \
                                                                          \
template<class T>                                                         \
inline MT::Array<T>& name##S(MT::Array<T>& x, T y, const MT::Array<T>& z){   \
  x.resizeAs(z);							  \
  T *xp=x.p;								  \
  const T *zp=z.p;							  \
  for(; xp!=x.pstop; xp++,zp++) *xp = y op *zp;                           \
  return x;                                                               \
}

BinaryOperation( plusA , + )
BinaryOperation( minusA , - )
BinaryOperation( mult , * )
BinaryOperation( div , / )
#undef BinaryOperation


//---------- compound assignment operators



//---------- unary functions

#define UnaryFunction( func )                                             \
template<class T>                                                         \
inline MT::Array<T> func (const MT::Array<T>& y){                         \
  MT::Array<T> x;                                                         \
  x.resizeAs(y);                                                          \
  T *xp=x.p;								  \
  const T *yp=y.p;							  \
  for(; xp!=x.pstop; xp++,yp++) *xp = ::func( *yp );                      \
  x.temp=true; return x;                                                  \
}									  \

// trigonometric functions
UnaryFunction( acos )
UnaryFunction( asin )
UnaryFunction( atan )
UnaryFunction( cos )
UnaryFunction( sin )
UnaryFunction( tan )

// hyperbolic functions
UnaryFunction( cosh )
UnaryFunction( sinh )
UnaryFunction( tanh )
UnaryFunction( acosh )
UnaryFunction( asinh )
UnaryFunction( atanh )

// exponential and logarithmic functions
UnaryFunction( exp )
UnaryFunction( log )
UnaryFunction( log10 )

// power functions
UnaryFunction( sqrt )
UnaryFunction( cbrt )

// nearest integer and absolute value
UnaryFunction( ceil )
UnaryFunction( fabs )
UnaryFunction( floor )
#undef UnaryFunction


//---------- binary functions

#define BinaryFunction( func )                                            \
template<class T>                                                         \
inline MT::Array<T> func(const MT::Array<T>& y,const MT::Array<T>& z){    \
  CHECK(y.N==z.N,							  \
  "binary operator on different array dimensions ("<<y.N<<", "<<z.N<<")");\
  MT::Array<T> x;                                                         \
  x.resizeAs(y);                                                            \
  for(uint i=x.N;i--; ) x.p[i]= func(y.p[i],z.p[i]);			  \
  x.temp=true; return x;                                                  \
}                                                                         \
                                                                          \
template<class T>                                                         \
inline MT::Array<T> func(const MT::Array<T>& y,T z){                      \
  MT::Array<T> x;                                                         \
  x.resizeAs(y);                                                            \
  for(uint i=x.N;i--; ) x.p[i]= func(y.p[i],z);			          \
  x.temp=true; return x;                                                  \
}                                                                         \
                                                                          \
template<class T>                                                         \
inline MT::Array<T> func(T y,const MT::Array<T>& z){                      \
  MT::Array<T> x;                                                         \
  x.resizeAs(z);                                                            \
  for(uint i=x.N;i--; ) x.p[i]= func(y,z.p[i]);			          \
  x.temp=true; return x;                                                  \
}

BinaryFunction( atan2 )
BinaryFunction( pow )
BinaryFunction( fmod )
#undef BinaryFunction

#ifndef MT_doxy // exclude these macros when generating the documentation
                // (doxygen can't handle them...)
#endif //(doxygen exclusion)


     /*

//! element-wise linear combination (plus with scalar factors for each array)
template<class T>
void plusSASA(MT::Array<T>& x,T a,const MT::Array<T>& y,T b,const MT::Array<T>& z){
  CHECK(y.N==z.N,"must have same size for adding!");
  uint i,n=y.N;
  x.resizeAs(y);
  for(i=0;i<n;i++) x.p[i]=a*y.p[i]+b*z.p[i];
  }*/


//===========================================================================
//
//!@name matrix operations
//

//! returns the identity matrix
inline MT::Array<double> Identity(uint dim){
  MT::Array<double> z;
  z.setId(dim);
  return z;
}

//! make symmetric \f$A=(A+A^T)/2\f$
inline void makeSymmetric(MT::Array<double>& A){
  CHECK(A.nd==2 && A.d0==A.d1,"not symmetric");
  uint n=A.d0,i,j;
  for(i=1;i<n;i++) for(j=0;j<i;j++) A(j,i) = A(i,j) = .5 * (A(i,j) + A(j,i));
}

//! make its transpose \f$A \gets A^T\f$
inline void transpose(MT::Array<double>& A){
  CHECK(A.nd==2 && A.d0==A.d1,"not symmetric");
  uint n=A.d0,i,j;
  double z;
  for(i=1;i<n;i++) for(j=0;j<i;j++){ z=A(j,i); A(j,i)=A(i,j); A(i,j)=z; }
}

//! x = y^T
template<class T>
void transpose(MT::Array<T>& x,const MT::Array<T>& y){
  CHECK(y.nd<=3,"can only transpose up to 3D arrays");
  if(y.nd==3){
    uint i,j,k,d0=y.d2,d1=y.d1,d2=y.d0;
    x.resize(d0,d1,d2);
    for(i=0;i<d0;i++) for(j=0;j<d1;j++) for(k=0;k<d2;k++)
      x(i,j,k) = y(k,j,i);
      //x.p[(i*d1+j)*d2+k]=y.p[(k*d1+j)*d0+i];
    return;
  }
  if(y.nd==2){
    uint i,j,d0=y.d1,d1=y.d0;
    x.resize(d0,d1);
    for(i=0;i<d0;i++) for(j=0;j<d1;j++) x.p[i*d1+j]=y.p[j*d0+i];
    return;
  }
  if(y.nd==1){
    x=y;
    x.reshape(1,y.N);
    return;
  }
  HALT("transpose not implemented for this dims");
}

//! returns the diagonal x = diag(y) (the diagonal-vector of the symmetric 2D matrix y)
template<class T>
inline void diag(MT::Array<T>& x,const MT::Array<T>& y){
  CHECK(y.nd==2 && y.d0==y.d1,"can only give diagonal of symmetric 2D matrix");
  x.resize(y.d0);
  uint i;
  for(i=0;i<x.d0;i++) x(i)=y(i,i);
}

//! sets x to be the diagonal matrix with diagonal v
template<class T>
inline void setDiagonal(MT::Array<T>& x,const MT::Array<T>& v){
  CHECK(v.nd==1,"can only give diagonal of 1D array");
  x.resize(v.d0,v.d0);
  x.setZero();
  uint i;
  for(i=0;i<v.d0;i++) x(i,i)=v(i);
}

//! inverse of a 2d matrix
template<class T>
void inverse2d(MT::Array<T>& invA,const MT::Array<T>& A){
  invA.resize(2,2);
  invA(0,0)=A(1,1); invA(1,1)=A(0,0); invA(0,1)=-A(0,1); invA(1,0)=-A(1,0);
  invA/=A(0,0)*A(1,1)-A(0,1)*A(1,0);
}

//! constructs the block matrix X=[A,B ; C,D]
template<class T>
void blockMatrix(MT::Array<T>& X,const MT::Array<T>& A,const MT::Array<T>& B,const MT::Array<T>& C,const MT::Array<T>& D){
  CHECK(A.nd==2 && B.nd==2 && C.nd==2 && D.nd==2,"");
  CHECK(A.d0==B.d0 && A.d1==C.d1 && B.d1==D.d1 && C.d0==D.d0,"");
  uint i,j,a=A.d0,b=A.d1;
  X.resize(A.d0+C.d0,A.d1+B.d1);
  for(i=0;i<A.d0;i++) for(j=0;j<A.d1;j++) X(i  ,j  )=A(i,j);
  for(i=0;i<B.d0;i++) for(j=0;j<B.d1;j++) X(i  ,j+b)=B(i,j);
  for(i=0;i<C.d0;i++) for(j=0;j<C.d1;j++) X(i+a,j  )=C(i,j);
  for(i=0;i<D.d0;i++) for(j=0;j<D.d1;j++) X(i+a,j+b)=D(i,j);
}

//! constructs a vector x=[a,b]
template<class T>
void blockVector(MT::Array<T>& x,const MT::Array<T>& a,const MT::Array<T>& b){
  CHECK(a.nd==1 && b.nd==1,"");
  uint i;
  x.resize(a.N+b.N);
  for(i=0;i<a.N;i++) x(i    )=a(i);
  for(i=0;i<b.N;i++) x(i+a.N)=b(i);
}

//! constructs the block matrix X=[A,B ; C,D] where D is a scalar
template<class T>
void blockMatrix(MT::Array<T>& X,const MT::Array<T>& A,const MT::Array<T>& B,const MT::Array<T>& C,const double& D){
  CHECK(A.nd==2 && B.nd==1 && C.nd==1,"");
  CHECK(A.d0==B.d0 && A.d1==C.d0,"");
  uint i,j,a=A.d0,b=A.d1;
  X.resize(A.d0+1,A.d1+1);
  for(i=0;i<a;i++) for(j=0;j<b;j++) X(i,j)=A(i,j);
  for(i=0;i<a;i++) X(i,b)=B(i);
  for(j=0;j<b;j++) X(a,j)=C(j);
  X(a,b)=D;
}


//===========================================================================
//
//!@name tensor operations
//

#define maxRank 10

/*! \f$Y_{i_Yid(0),i_Yid(1)} = \sum_{i_1} X_{i_0,i_1,i_2}\f$. Get the marginal Y
  from X, where Y will share the slots `Yid' with X */
template<class T>
void tensorMarginal2(MT::Array<T> &Y,const MT::Array<T> &X,const uintA &Yid){
  uint i,j,k;
  int l;
  CHECK(Yid.N<=X.nd,"can't take slots "<<Yid<<" from "<<X.nd<<"D tensor");
  CHECK(Yid.N>0    ,"can't take 0 slots (tensorMarginal2) from "<<X.nd<<"D tensor");

  //loop dimensions
  uint N=X.N,rankX=X.nd;
  uint D[maxRank];  for(i=0;i<rankX;i++) D[i]=X.d[i];
  uint I[maxRank];  for(i=0;i<rankX;i++) I[i]=0;  //index tuple

  //Y dimensions, permutation & increments
  uint rankY=Yid.N;
  uint Ydim[maxRank];   for(i=0;i<rankY;i++) Ydim[i]=D[Yid.p[i]];
  Y.resize(rankY,Ydim);
  int  Yslot[maxRank];  for(i=0;i<rankX;i++) Yslot[i]=-1;  for(i=0;i<rankY;i++) Yslot[Yid.p[i]]=i;   
  uint Yinc[maxRank];   Yinc[rankY-1]=1;  for(i=rankY-1;i--;) Yinc[i] = Ydim[i+1] * Yinc[i+1];
  
  Y.setZero();
  for(i=0,j=0;;){
    Y.p[j] += X.p[i];
    //increment all indices
    i++;
    if(i==N) break; //loop exit
    for(k=rankX;k--;){
      I[k]++;
      l=Yslot[k]; if(l!=-1) j+=Yinc[l]; //increment Y index
      if(I[k]<D[k]) break;
      I[k]=0;
      if(l!=-1) if(!l) j-=Y.N; else j-=Yinc[l-1];
    }
#if 0
    //check if incrementing Y worked out
    uint jj=0;
    for(k=0;k<Yid.N;k++){ jj*=Ydim(k); jj+=I(Yid(k)); }
    CHECK(jj==j,"");
    //check if incrementing I worked out
    uintA II;
    getIndexTuple(II,i,D);
    CHECK(II==I,"not equal: "<<II <<I);
#endif
  }
}

/*! \f$Y_{i_2,i_0} = {\rm max}_{i_1} X_{i_0,i_1,i_2}\f$. Get the ``max-marginal'' Y
  from X, where Y will share the slots `Yid' with X (basis of max-product BP) */
template<class T>
void tensorMaxMarginal2(MT::Array<T> &Y,const MT::Array<T> &X,const uintA &Yid){
  uint i,j,k;
  int l;
  CHECK(Yid.N<=X.nd,"can't take slots "<<Yid<<" from "<<X.nd<<"D tensor");

  //loop dimensions
  uint N=X.N,rankX=X.nd;
  uint D[maxRank];  for(i=0;i<rankX;i++) D[i]=X.d[i];
  uint I[maxRank];  for(i=0;i<rankX;i++) I[i]=0;  //index tuple

  //Y dimensions, permutation & increments
  uint rankY=Yid.N;
  uint Ydim[maxRank];   for(i=0;i<rankY;i++) Ydim[i]=D[Yid.p[i]];
  Y.resize(rankY,Ydim);
  int  Yslot[maxRank];  for(i=0;i<rankX;i++) Yslot[i]=-1;  for(i=0;i<rankY;i++) Yslot[Yid.p[i]]=i;   
  uint Yinc[maxRank];   Yinc[rankY-1]=1;  for(i=rankY-1;i--;) Yinc[i] = Ydim[i+1] * Yinc[i+1];

  Y.setZero();
  for(i=0,j=0;;){
    if(X.p[i]>Y.p[j]) Y.p[j] = X.p[i];
    //increment all indices
    i++;
    if(i==N) break; //loop exit
    for(k=rankX;k--;){
      I[k]++;
      l=Yslot[k]; if(l!=-1) j+=Yinc[l]; //increment Y index
      if(I[k]<D[k]) break;
      I[k]=0;
      if(l!=-1) if(!l) j-=Y.N; else j-=Yinc[l-1];
    }
#if 0
    //check if incrementing Y worked out
    uint jj=0;
    for(k=0;k<rankY;k++){ jj*=Ydim(k); jj+=I(Yid(k)); }
    CHECK(jj==j,"");
    //check if incrementing I worked out
    uintA II;
    getIndexTuple(II,i,D);
    CHECK(II==I,"not equal: "<<II <<I);
#endif
  }
}

template<class T>
void tensorMarginal_old(MT::Array<T> &y,const MT::Array<T> &x,const uintA &xd,const uintA &ids){
  uint i,j,k,n=product(xd);
  CHECK(x.N==n,"");
  //size y
  uintA yd(ids.N);
  for(i=0;i<ids.N;i++) yd(i)=xd(ids(i));
  //y.tensorResize(yd); y.setZero();
  y.resize(product(yd)); y.setZero();

  uintA xt(xd.N); xt.setZero();
  for(i=0;i<n;i++){
    //compute j
    for(j=0,k=0;k<ids.N;k++){ j*=yd(k); j+=xt(ids(k)); }
    //uintA yt(ids.N); for(k=ids.N;k--;){ yt(k)=xt(ids(k)); }
    //cout <<"i=" <<i <<" j=" <<j <<" yt=" <<yt <<" xt=" <<xt <<endl;
    y(j)+=x.elem(i);
    //increment xt
    for(k=xt.N;k--;){ xt(k)++; if(xt(k)<xd(k)) break; else xt(k)=0; }
  }
}

/*! \f$X_{i_0,i_1,i_2} \gets X_{i_0,i_1,i_2}~ Y_{i_Yid(0),i_Yid(1)}\f$. Multiply Y onto X,
  where Y shares the slots `Yid' with X */
template<class T>
void tensorMultiply2(MT::Array<T> &X,const MT::Array<T> &Y,const uintA &Yid){
  uint i,j,k;
  int l;
  CHECK(Yid.N==Y.nd,"need to specify "<<Y.nd<<" slots, not "<<Yid.N);
  CHECK(Y.nd<=X.nd,"can't multiply "<<Y.nd<<" slots onto "<<X.nd<<"D tensor");

  //loop dimensions
  uint N=X.N,rankX=X.nd;
  uint D[maxRank];  for(i=0;i<rankX;i++) D[i]=X.d[i];
  uint I[maxRank];  for(i=0;i<rankX;i++) I[i]=0;  //index tuple

  //Y dimensions, permutation & increments
  uint rankY=Yid.N;
  uint Ydim[maxRank];  for(i=0;i<rankY;i++) Ydim[i]=D[Yid.p[i]];
  CHECK(uintA(Ydim,rankY)==Y.d,"tensorMultiply: mismatch in dimensions: Xdim=" <<D <<"Ydim=" <<uintA(Y.d,Y.nd) <<"Yid="<<Yid <<"Xdim(Yid)="<<uintA(Ydim,rankY));
  int  Yslot[maxRank];    for(i=0;i<rankX;i++) Yslot[i]=-1;  for(i=0;i<rankY;i++) Yslot[Yid.p[i]]=i;   
  uint Yinc[maxRank];  Yinc[rankY-1]=1;  for(i=rankY-1;i--;) Yinc[i] = Ydim[i+1] * Yinc[i+1];

  for(i=0,j=0;;){
    X.p[i] *= Y.p[j];
    //increment all indices
    i++;
    if(i==N) break; //loop exit
    for(k=rankX;k--;){
      I[k]++;
      l=Yslot[k]; if(l!=-1) j+=Yinc[l]; //increment Y index
      if(I[k]<D[k]) break;
      I[k]=0;
      if(l!=-1) if(!l) j-=Y.N; else j-=Yinc[l-1];
    }
#if 0
    //check if incrementing Y worked out
    uint jj=0;
    for(k=0;k<rankY;k++){ jj*=Ydim(k); jj+=I(Yid(k)); }
    CHECK(jj==j,"");
    //check if incrementing I worked out
    uintA II;
    getIndexTuple(II,i,D);
    CHECK(II==I,"not equal: "<<II <<I);
#endif
  }
}

/*! multiply y onto x, where x has dimensions `d' and y shares the
  dimensions `ids' with x */
template<class T>
void tensorMultiply_old(MT::Array<T> &x,const MT::Array<T> &y,const uintA &d,const uintA &ids){
  uint i,j,k,n=x.N;
  CHECK(n==product(d),"");

  uintA yd(ids.N);
  for(i=0;i<ids.N;i++) yd(i)=d(ids(i));
  CHECK(y.N==product(yd),"");

  uintA I(d.N); I.setZero();
  for(i=0;i<n;i++){
    for(j=0,k=0;k<ids.N;k++){ j*=yd(k); j+=I(ids(k)); }
    x.elem(i) *= y.elem(j);
    for(k=I.N;k--;){ I(k)++; if(I(k)<d(k)) break; else I(k)=0; }
  }
}


inline uintA TUP(uint i){                                           uintA z(1); z(0)=i; return z; }
inline uintA TUP(uint i,uint j){                                    uintA z(2); z(0)=i; z(1)=j; return z; }
inline uintA TUP(uint i,uint j,uint k){                             uintA z(3); z(0)=i; z(1)=j; z(2)=k; return z; }
inline uintA TUP(uint i,uint j,uint k,uint l){                      uintA z(4); z(0)=i; z(1)=j; z(2)=k; z(3)=l; return z; }
inline uintA TUP(uint i,uint j,uint k,uint l,uint m){               uintA z(5); z(0)=i; z(1)=j; z(2)=k; z(3)=l; z(4)=m; return z; }
inline uintA TUP(uint i,uint j,uint k,uint l,uint m,uint n){        uintA z(6); z(0)=i; z(1)=j; z(2)=k; z(3)=l; z(4)=m; z(5)=n; return z; }
inline uintA TUP(uint i,uint j,uint k,uint l,uint m,uint n,uint o){ uintA z(7); z(0)=i; z(1)=j; z(2)=k; z(3)=l; z(4)=m; z(5)=n; z(6)=o; return z; }
inline uintA TUP(uint i,uint j,uint k,uint l,uint m,uint n,uint o,uint p){ uintA z(8); z(0)=i; z(1)=j; z(2)=k; z(3)=l; z(4)=m; z(5)=n; z(6)=o; z(7)=p; return z; }

inline arr ARR(double i){                                           arr z(1); z(0)=i; return z; }
inline arr ARR(double i,double j){                                    arr z(2); z(0)=i; z(1)=j; return z; }
inline arr ARR(double i,double j,double k){                             arr z(3); z(0)=i; z(1)=j; z(2)=k; return z; }
inline arr ARR(double i,double j,double k,double l){                      arr z(4); z(0)=i; z(1)=j; z(2)=k; z(3)=l; return z; }
inline arr ARR(double i,double j,double k,double l,double m){               arr z(5); z(0)=i; z(1)=j; z(2)=k; z(3)=l; z(4)=m; return z; }
inline arr ARR(double i,double j,double k,double l,double m,double n){        arr z(6); z(0)=i; z(1)=j; z(2)=k; z(3)=l; z(4)=m; z(5)=n; return z; }
inline arr ARR(double i,double j,double k,double l,double m,double n,double o){ arr z(7); z(0)=i; z(1)=j; z(2)=k; z(3)=l; z(4)=m; z(5)=n; z(6)=o; return z; }

void getIndexTuple(uintA &I,uint i,const uintA &d);

/*! X becomes a tensor product (maybe with some variables summed out)
    of A and B. Aid and Bid indicate which slots of A and B are
    associated with which slots of C and the summation sign. More
    precisely, if we have \f$C_{i_0i_1i_2} = \sum_{i_3i_4}
    A_{i_4i_2i_1}~ B_{i_3i_0}\f$ then you should call
    tensor(C,A,TUP(4,2,1),B,TUP(3,0),2); Here, the `2` indicates that
    the last two indices of i_0,..,i_4 are summed over, and C only
    becomes a 3rd rank instead of 5th rank tensor */
void tensor(arr &X,
	    const arr &A, const uintA &Aid,
	    const arr &B, const uintA &Bid,
	    uint sum=0);


//===========================================================================
//
//!@name probability distribution operations
//

/*!\brief entropy \f$H_i = - \sum_x p(X_i=x) \ln p(X_i=x)\f$, where \f$x\f$ ranges
  from \c 0 to \c range-1, of the \c ith variable */
template<class T>
inline T entropy(const MT::Array<T>& v){
  T t(0);
  for(uint i=v.N; i--;) if(v.p[i]) t-=v.p[i]*::log(v.p[i]);
  return t/MT_LN2;
}

//! v = v / sum(v)
template<class T>
inline T normalizeDist(MT::Array<T>& v){
  double Z=sum(v);
  if(Z>1e-100) v/=Z; else v=1./(double)v.N;
  return Z;
}

//! v = v / sum(v)
template<class T>
inline void makeConditional(MT::Array<T>& P){
  MT_MSG("makeConditional: don't use this anymore because it normalizes over the second index!!!, rather use tensorCondNormalize and condition on _later_ indices");
  CHECK(P.nd==2,"");
  uint i,j;
  double pi;
  for(i=0;i<P.d0;i++){
    pi=0.;
    for(j=0;j<P.d1;j++) pi+=P(i,j);
    for(j=0;j<P.d1;j++) P(i,j) /= pi;
  }
}

//! check whether this is a distribution over the first index w.r.t. the later indices
template<class T>
inline void checkNormalization(MT::Array<T>& v,double tol=1e-10){
  double p;
  uint i,j,k;
  switch(v.nd){
  case 1:
    for(p=0,i=0;i<v.d0;i++) p+=v(i);
    CHECK(fabs(1.-p)<tol,"distribution is not normalized: "<<v);
    break;
  case 2:
    for(j=0;j<v.d1;j++){
      for(p=0,i=0;i<v.d0;i++) p+=v(i,j);
      CHECK(fabs(1.-p)<tol,"distribution is not normalized: "<<v);
    }
    break;
  case 3:
    for(j=0;j<v.d1;j++) for(k=0;k<v.d2;k++){
      for(p=0,i=0;i<v.d0;i++) p+=v(i,j,k);
      CHECK(fabs(1.-p)<tol,"distribution is not normalized: "<<v);
    }
    break;
  }
}

template<class T>
inline void eliminate(MT::Array<T>& x,const MT::Array<T>& y,uint d){
  CHECK(y.nd==2,"only implemented for 2D yet");
  uint i,j;
  if(d==1){
    x.resize(y.d0); x=0.;
    for(i=0;i<y.d0;i++) for(j=0;j<y.d1;j++) x(i)+=y(i,j);
  }
  if(d==0){
    x.resize(y.d1); x=0.;
    for(i=0;i<y.d0;i++) for(j=0;j<y.d1;j++) x(j)+=y(i,j);
  }
}

template<class T>
inline void eliminate(MT::Array<T>& x,const MT::Array<T>& y,uint d,uint e){
  CHECK(y.nd==3,"only implemented for 3D yet");
  uint i,j,k;
  if(d==1 && e==2){
    x.resize(y.d0); x=0.;
    for(i=0;i<y.d0;i++) for(j=0;j<y.d1;j++) for(k=0;k<y.d2;k++) x(i)+=y(i,j,k);
  }
  if(d==0 && e==2){
    x.resize(y.d1); x=0.;
    for(i=0;i<y.d0;i++) for(j=0;j<y.d1;j++) for(k=0;k<y.d2;k++) x(j)+=y(i,j,k);
  }
  if(d==0 && e==1){
    x.resize(y.d2); x=0.;
    for(i=0;i<y.d0;i++) for(j=0;j<y.d1;j++) for(k=0;k<y.d2;k++) x(k)+=y(i,j,k);
  }
}

// Eliminates one-dimension, d, from a 3D-tensor, y, and puts the result in x.
template<class T>
inline void eliminatePartial(MT::Array<T>& x,const MT::Array<T>& y,uint d){
  CHECK(y.nd==3,"only implemented for 3D yet");
  uint i,j,k;
  if(d==2){
    x.resize(y.d0,y.d1); x=0.;
    for(i=0;i<y.d0;i++) for(j=0;j<y.d1;j++) for(k=0;k<y.d2;k++) x(i,j)+=y(i,j,k);
  }
  if(d==1){
    x.resize(y.d0,y.d2); x=0.;
    for(i=0;i<y.d0;i++) for(j=0;j<y.d1;j++) for(k=0;k<y.d2;k++) x(i,k)+=y(i,j,k);
  }
  if(d==0){
    x.resize(y.d1,y.d2); x=0.;
    for(i=0;i<y.d0;i++) for(j=0;j<y.d1;j++) for(k=0;k<y.d2;k++) x(j,k)+=y(i,j,k);
  }
}

/*! stochastic universal sampling: given a probability vector p
  (that is positive in each component and sums to one), the returned
  array s is resized to length n and contains the states of n
  samples. */
void SUS(const arr& p,uint n,uintA& s);

/*! as above, but returns only the index of one sample */
uint SUS(const arr& p);


//===========================================================================
//
//!@name distances
//

//! \f$\sum_i (v^i-w^i)^2\f$
template<class T> 
inline T sqrDistance(const MT::Array<T>& v, const MT::Array<T>& w){
  CHECK(v.N==w.N, 
    "sqrDistance on different array dimensions ("<<v.N<<", "<<w.N<<")");
  T d,t(0);
  for(uint i=v.N;i--;){ d=v.p[i]-w.p[i]; t+=d*d; }
  return t;
}

template<class T>
inline T maxDiff(const MT::Array<T>& v,const MT::Array<T>& w,uint *im=0){
  CHECK(v.N==w.N, 
    "maxDiff on different array dimensions ("<<v.N<<", "<<w.N<<")");
  T d,t(0);
  if(!im) for(uint i=v.N;i--;){ d=fabs(v.p[i]-w.p[i]); if(d>t) t=d; }
  else    for(uint i=v.N;i--;){ d=fabs(v.p[i]-w.p[i]); if(d>t){ t=d; *im=i; } }
  return t;
}

template<class T>
inline T maxRelDiff(const MT::Array<T>& v,const MT::Array<T>& w,double tol){
  CHECK(v.N==w.N, 
    "maxDiff on different array dimensions ("<<v.N<<", "<<w.N<<")");
  T d,t(0),a,b,c;
  for(uint i=v.N;i--;){
    a=fabs(v.p[i]) + tol;
    b=fabs(w.p[i]) + tol;
    if(a<b){ c=a; a=b; b=c; }
    d=a/b-1.;
    if(d>t) t=d;
  }
  return t;
}

//! \f$\sum_{i|{\rm mask}_i={\rm true}} (v^i-w^i)^2\f$
template<class T> 
inline T sqrDistance(const MT::Array<T>& v, const MT::Array<T>& w,const MT::Array<bool>& mask){
  CHECK(v.N==w.N, 
    "sqrDistance on different array dimensions ("<<v.N<<", "<<w.N<<")");
  T d,t(0);
  for(uint i=v.N;i--;) if(mask(i)){ d=v.p[i]-w.p[i]; t+=d*d; }
  return t;
}

//! \f$\sqrt{\sum_{ij} g_{ij} (v^i-w^i) (v^j-w^j)}\f$
template<class T>
inline T sqrDistance(const MT::Array<T>& g,const MT::Array<T>& v,const MT::Array<T>& w){
  MT::Array<T> d;
  ::minusA(d,v,w);
  return scalarProduct(g,d,d);
}

//! \f$\sqrt{\sum_i (v^i-w^i)^2}\f$
template<class T> 
inline T euclideanDistance(const MT::Array<T>& v, const MT::Array<T>& w){
  return ::sqrt(sqrDistance(v,w));
}

//! \f$\sqrt{\sum_i (v^i-w^i)^2}\f$
template<class T> 
inline T metricDistance(const MT::Array<T>& g,const MT::Array<T>& v, const MT::Array<T>& w){
  return ::sqrt(sqrDistance(g,v,w));
}


//===========================================================================
//
//!@name running sums
//

//! \f$\sum_i x_i\f$
template<class T>
inline T sum(const MT::Array<T>& v){
  T t(0);
  for(uint i=v.N; i--; t+=v.p[i]);
  return t;
}

//! \f$\sum_i |x_i|\f$
template<class T>
inline T sumOfAbs(const MT::Array<T>& v){
  T t(0);
  for(uint i=v.N; i--; t+=fabs(v.p[i]));
  return t;
}

//! \f$\sum_i x_i^2\f$
template<class T>
inline T sumOfSqr(const MT::Array<T>& v){
  T t(0);
  for(uint i=v.N; i--; t+=v.p[i]*v.p[i]);
  return t;
}

//! \f$\sqrt{\sum_i x_i^2}\f$
template<class T>
inline T norm(const MT::Array<T>& v){ return sqrt(sumOfSqr(v)); }

//! \f$\sqrt{\sum_i x_i^2}\f$
template<class T>
inline T mean(const MT::Array<T>& v){ return sum(v)/v.N; }

template<class T>
inline T var(const MT::Array<T>& v){ double m=mean(v); return (sumOfSqr(v)-m*m/v.N)/v.N; }

//! \f$\sum_i x_{ii}\f$
template<class T>
inline T trace(const MT::Array<T>& v){
  CHECK(v.nd==2 && v.d0==v.d1,"only for squared matrix");
  T t(0);
  for(uint i=0;i<v.d0;i++) t+=v(i,i);
  return t;
}

template<class T>
inline T minDiag(const MT::Array<T>& v){
  CHECK(v.nd==2 && v.d0==v.d1,"only for squared matrix");
  T t=v(0,0);
  for(uint i=1;i<v.d0;i++) if(v(i,i)<t) t=v(i,i);
  return t;
}


//===========================================================================
//
//!@name cblas and Lapack support
//

namespace MT{
//! use this to turn on Lapack routines [default true if MT_LAPACK is defined]
extern bool useLapack;
}

#ifdef MT_LAPACK
extern "C"{
#include"cblas.h"
#include"f2c.h"
#include"clapack.h"
}

//matrix-matrix multiplication
inline void blasMM(void *pX,void *pA,void *pB){
  arr &X=*(arr*)pX,&A=*(arr*)pA,&B=*(arr*)pB;
  CHECK(A.d1==B.d0,"matrix multiplication: wrong dimensions");
  X.resize(A.d0,B.d1);
  cblas_dgemm(CblasRowMajor,
	      CblasNoTrans,CblasNoTrans,
	      A.d0,B.d1,A.d1,
	      1.,A.p,A.d1,
	      B.p,B.d1,1.,
	      X.p,X.d1);
#ifdef MT_Linux
  //bad joke? the 10th element is always nan...
  if(X.N>9){
    X(0,9)=0.; for(uint k=0;k<A.d1;k++) X(0,9) += A(0,k) * B(k,9);
  }
#endif

#if 0 //test
  arr Y(A.d0,B.d1);
  uint i,j,k;
  Y.setZero();
  for(i=0;i<Y.d0;i++) for(j=0;j<Y.d1;j++) for(k=0;k<A.d1;k++)
    Y(i,j) += A(i,k) * B(k,j);
  std::cout  <<"blasMM error = " <<sqrDistance(X,Y) <<std::endl;
#endif
}

inline uint lapackSVD(const arr& A,
		      arr& U,
		      arr& d,
		      arr& Vt){
  static arr Atmp,work;
  Atmp=A;
  //transpose(Atmp,A);
  integer M=A.d0,N=A.d1,D=M<N?M:N;
  U.resize(M,D);
  d.resize(D);
  Vt.resize(D,N);
  work.resize(10*(M+N));
  integer info,wn=work.N;
  dgesvd_("S", "S", &N, &M, Atmp.p, &N,
	  d.p, Vt.p, &N, U.p, &D,
	  work.p,&wn, &info);
  if(info){
    std::cerr <<"LAPACK SVD error info = " <<info <<std::endl;
  }

#if 0 //test!
  MT::useLapack=false;
  arr dD,I;
  setDiagonal(dD,d);
  //cout <<U <<dD <<Vt;
  //Atmp = V * D * U;
  Atmp = U * dD * Vt;
  std::cout <<"SVD is correct:  " <<sqrDistance(Atmp,A) <<' ' <<endl;
  setIdentity(I,D);
  std::cout <<"U is orthogonal: " <<sqrDistance(~U * U,I) <<' ' <<endl;
  std::cout <<"Vt is orthogonal: " <<sqrDistance(Vt * ~Vt,I) <<endl;
  MT::useLapack=true;
#endif
  return D;
}
#endif


//===========================================================================
//
//!@name SVD etc
//

//! called from svd if MT_LAPACK is not defined
uint own_svd_implementation(const arr& A,
			    arr& U,
			    arr& w,
			    arr& V,
			    bool sort);

/*!\brief Singular Value Decomposition (from Numerical Recipes);
  computes \f$U, D, V\f$ with \f$A = U D V^T\f$ from \f$A\f$ such that
  \f$U\f$ and \f$V\f$ are orthogonal and \f$D\f$ diagonal (the
  returned array d is 1-dimensional) -- uses LAPACK if MT_LAPACK is
  defined */
inline uint svd(const arr& A,arr& U,arr& d,arr& V,bool sort=true){
#ifdef MT_LAPACK
  if(MT::useLapack) return lapackSVD(A,U,d,V);
  V=~V;
#endif
  uint r=own_svd_implementation(A,U,d,V,sort);
  return r;
}

//! gives a decomposition \f$A = U V^T\f$
inline void svd(const arr& A,arr& U,arr& V){
  arr d,D;
  ::svd(A,U,d,V);
  D.resize(d.N,d.N); D=0.;
  for(uint i=0;i<d.N;i++) D(i,i)=::sqrt(d(i));
  U=U*D;
  V=V*D;
  //CHECK(maxDiff(A,U*~V) <1e-4,"");
}

//! Pseudo Inverse based on SVD; computes \f$B\f$ such that \f$ABA = A\f$
uint inverse_SVD(arr& inverse,const arr& A);

namespace MT{  void inverse_LU(arr& invX,const arr& X); }

inline uint inverse(arr& inverse,const arr& A){
  return inverse_SVD(inverse,A);
  //MT::inverse_LU(inverse,A); return A.d0;
}

//! weighted pseudo inverse
void pseudoInverse(arr& inverse,const arr& A,const arr& W,double robustnessEps=0.);

//! calls inverse(B,A) and returns B
inline arr inverse(const arr& A){ arr B; inverse(B,A); return B; }

//! the determinant of a 2D squared matrix
double determinant(const arr& A);

/*!\brief the cofactor is the determinant of a 2D squared matrix after removing
  the ith row and the jth column */
double cofactor(const arr& A,uint i,uint j);


//===========================================================================
//
//!@name products
//

//! \f$\prod_i x_i\f$
template<class T>
inline T product(const MT::Array<T>& v){
  T t(1);
  for(uint i=v.N; i--; t *= v.p[i]);
  return t;
}

/*!\brief inner product (also ordinary matrix or scalar product):
  \f$\forall_{ik}:~ x_{ik} = \sum_j v_{ij}\, w_{jk}\f$ but also:
  \f$\forall_{i}:~ x_{i} = \sum_j v_{ij}\, w_{j}\f$*/
template<class T> 
void innerProduct(MT::Array<T>& x,const MT::Array<T>& y, const MT::Array<T>& z){
  /*
  if(y.nd==2 && z.nd==2 && y.N==z.N && y.d1==1 && z.d1==1){  //elem-wise
    HALT("make element-wise multiplication explicite!");
    mult(x,y,z);
    return;
  }
  */
  if(y.nd==2 && z.nd==3){
    MT::Array<T> zz; zz.referTo(z);
    zz.reshape(z.d0,z.d1*z.d2);
    innerProduct(x,y,zz);
    x.reshape(y.d0,z.d1,z.d2);
    return;
  }
  if(y.nd==3 && z.nd==2){
    MT::Array<T> yy; yy.referTo(y);
    yy.reshape(y.d0*y.d1,y.d2);
    innerProduct(x,yy,z);
    x.reshape(y.d0,y.d1,z.d1);
    return;
  }
  if(y.nd==3 && z.nd==1){
    MT::Array<T> yy; yy.referTo(y);
    yy.reshape(y.d0*y.d1,y.d2);
    innerProduct(x,yy,z);
    x.reshape(y.d0,y.d1);
    return;
  }
  if(y.nd==1 && z.nd==3){
    MT::Array<T> zz; zz.referTo(z);
    zz.reshape(z.d0,z.d1*z.d2);
    innerProduct(x,y,zz);
    x.reshape(z.d1,z.d2);
    return;
  }
  if(y.nd==2 && z.nd==2){ //plain matrix multiplication
#ifdef MT_LAPACK
    if(MT::useLapack && typeid(T)==typeid(double)){ blasMM(&x,(void*)&y,(void*)&z); return; }
#endif
    CHECK(y.d1==z.d0,"wrong dimensions for inner product");
    uint i,j,d0=y.d0,d1=z.d1,dk=y.d1;
#if 0
    if(y.mtype==MT::Array<T>::diag){
      x.resize(d0,d1);
      for(i=0;i<d0;i++) for(j=0;j<d1;j++) x(i,j) = y(i,i) * z(i,j);
      return;
    }
    if(z.mtype==MT::Array<T>::diag){
      x.resize(d0,d1);
      for(i=0;i<d0;i++) for(j=0;j<d1;j++) x(i,j) *= y(i,j) * z(j,j);
      return;
    }
#endif
    T *a,*astop,*b,*c;
    x.resize(d0,d1); x.setZero();
    c=x.p;
    for(i=0;i<d0;i++) for(j=0;j<d1;j++){
      //for(s=0.,k=0;k<dk;k++) s+=y.p[i*dk+k]*z.p[k*d1+j];
      //this is faster:
      a=y.p+i*dk; astop=a+dk; b=z.p+j;
      for(;a!=astop; a++, b+=d1) (*c)+=(*a) * (*b);
      c++;
    }
    return;
  }
  if(y.nd==2 && z.nd==1){ //matrix x vector -> vector
    CHECK(y.d1==z.d0,"wrong dimensions for inner product");
    uint i,d0=y.d0,dk=y.d1;
    T *a,*astop,*b,*c;
    x.resize(d0); x.setZero();
    c=x.p;
    for(i=0;i<d0;i++){
      //for(s=0.,k=0;k<dk;k++) s+=y.p[i*dk+k]*z.p[k];
      //this is faster:
      a=y.p+i*dk; astop=a+dk; b=z.p;
      for(;a!=astop; a++, b++) (*c)+=(*a) * (*b);
      c++;
    }
    return;
  }
  if(y.nd==1 && z.nd==2 && z.d0==1){ //vector x vector^T -> matrix (outer product)
    CHECK(y.d0==z.d1,"wrong dimensions for inner product");
    uint i,j,d0=y.d0,d1=z.d1;
    x.resize(d0,d1);
    for(i=0;i<d0;i++) for(j=0;j<d1;j++) x(i,j)=y(i)*z(0,j);
    return;
  }
  if(y.nd==1 && z.nd==2){ //vector^T x matrix -> vector^T
    CHECK(y.d0==z.d0,"wrong dimensions for inner product");
    uint i,k,d0=z.d1,dk=y.d0;
    x.resize(d0);
    T s;
    for(i=0;i<d0;i++){
      for(s=0,k=0;k<dk;k++) s+=y.p[k]*z.p[k*d0+i];
      x.p[i]=s;
    }
    return;
  }
  if(y.nd==1 && z.nd==1 && z.N==1){ //vector multiplied with scalar (disguised as 1D vector)
    uint k,dk=y.N;
    x.resize(y.N);
    for(k=0;k<dk;k++) x.p[k]=y.p[k]*z.p[0];
    return;
  }
  if(y.nd==1 && z.nd==1){ //should be scalar product, but be careful
    HALT("what do you want? scalar product or element wise multiplication?");
    CHECK(y.d0==z.d0,"wrong dimensions for inner product");
    uint k,dk=y.d0;
    x.resize(1);
    T s;
    for(s=0,k=0;k<dk;k++) s+=y.p[k]*z.p[k];
    x.p[0]=s;
    return;
  }
  HALT("inner product - not yet implemented for these dimensions");
}

/*!\brief outer product (also exterior or tensor product): \f$\forall_{ijk}:~
  x_{ijk} = v_{ij}\, w_{k}\f$ */
template<class T> 
void outerProduct(MT::Array<T>& x,const MT::Array<T>& y, const MT::Array<T>& z){
  if(y.nd==1 && z.nd==1){
    uint i,j,d0=y.d0,d1=z.d0;
    x.resize(d0,d1);
    for(i=0;i<d0;i++) for(j=0;j<d1;j++) x.p[i*d1+j]=y.p[i]*z.p[j];
    return;
  }
  HALT("outer product - not yet implemented for these dimensions");
}

//! \f$\sum_i v_i\, w_i\f$, or \f$\sum_{ij} v_{ij}\, w_{ij}\f$, etc.
template<class T> 
inline T scalarProduct(const MT::Array<T>& v, const MT::Array<T>& w){
  CHECK(v.N==w.N, 
    "scalar product on different array dimensions ("<<v.N<<", "<<w.N<<")");
  T t(0);
  for(uint i=v.N; i--; t+=v.p[i]*w.p[i]);
  return t;
}

//! \f$\sum_{ij} g_{ij}\, v_i\, w_i\f$
template<class T> 
inline T scalarProduct(const MT::Array<T>& g,const MT::Array<T>& v, const MT::Array<T>& w){
  CHECK(v.N==w.N && g.nd==2 && g.d0==v.N && g.d1==w.N,
    "scalar product on different array dimensions ("<<v.N<<", "<<w.N<<")");
  T t(0);
  for(uint i=0;i<g.d0;i++) for(uint j=0;j<g.d1;j++) t+=g(i,j)*v.p[i]*w.p[j];
  return t;
}


//===========================================================================
//
//!@name set operations
//

//! x becomes the section of y and z
template<class T> 
inline void section(MT::Array<T>& x,const MT::Array<T>& y, const MT::Array<T>& z){
  uint i,j;
  x.clear();
  for(i=0;i<y.N;i++){
    for(j=0;j<z.N;j++) if(y(i)==z(j)){ x.append(y(i)); break; }
  }
}


//===========================================================================
//
//!@name min and max
//

//! \f$\min_i x_i\f$
template<class T>
inline T minA(const MT::Array<T>& v){
  CHECK(v.N>0,"");
  T t(v.p[0]);
  for(uint i=1; i<v.N; ++i) if(v.p[i]<t) t=v.p[i];
  return t;
}

//! also returns the index (argmin) \c ind
template<class T>
inline T minA(const MT::Array<T>& v, uint & ind, uint start=0, uint end=0){
  CHECK(v.N>0,"");
  CHECK(v.N>start,"");	
  CHECK(v.N>=end,"");
  CHECK(end>=start,"");
  T t(v(start));
  ind=start;
  if (end==0)	end=v.N;
  for(uint i=start; i<end; ++i) if(v.p[i]<t){
    t  =v.p[i];
    ind=i;
  }
  return t;
}

//! \f$\max_i x_i\f$
template<class T>
inline T maxA(const MT::Array<T>& v){
  CHECK(v.N>0,"");
  T t(v.p[0]);
  for(uint i=1; i<v.N; ++i) if(v.p[i]>t) t=v.p[i];
  return t;
}

//! also returns the index (argmax) \c ind
template<class T>
inline T maxA(const MT::Array<T>& v, uint & ind, uint start=0, uint end=0){
  CHECK(v.N>0,"");
  CHECK(v.N>start,"");
  CHECK(v.N>=end,"");
  CHECK(end>=start,"");
  T t(v(start));
  ind=start;
  if(end==0){
    end=v.N;
  }
  for(uint i=start; i<end; ++i) if(v.p[i]>t){
    t  =v.p[i];
    ind=long(i);
  }
  return t;
}

//! gets the min and max
template<class T>
inline void minmaxA(const MT::Array<T>& v, T& minVal, T& maxVal){
  if(v.N>0){
    minVal=maxVal=v.p[0];
    for(uint i=1; i<v.N; ++i){
      if(v.p[i]<minVal) minVal=v.p[i];
      else if(v.p[i]>maxVal) maxVal=v.p[i];
    }
  }
}

//! get absolute maximum (using fabs)
template<class T>
inline T absMax(const MT::Array<T>& x){
  CHECK(x.N>0,"");
  T t(x.p[0]);
  for(uint i=1; i<x.N; ++i) if(fabs(x.p[i])>t) t=fabs(x.p[i]);
  return t;
}


//===========================================================================
//
//!@name randomizations
//

//! Assign all elements of \c a to a uniformly distributed discrete value in {low,..,hi}
template<class T>
void rndInt(MT::Array<T>& a,int low=0.,int high=1.){
  for(uint i=0;i<a.N;i++) a.p[i]=low+(int)rnd.num(1+high-low);
}

//! Assign all elements of \c a to a uniformly distributed continuous value in [low,hi]
template<class T>
void rndUni(MT::Array<T>& a,double low=0.,double high=1.,bool add=false){
  if(!add) for(uint i=0;i<a.N;i++) a.p[i]=rnd.uni(low,high);
  else     for(uint i=0;i<a.N;i++) a.p[i]+=rnd.uni(low,high);
}

//! Assign all elements of \c a to a Gaussian random variable
template<class T>
void rndGauss(MT::Array<T>& a,double stdDev,bool add=false){
  if(a.nd==2 && a.d0!=a.d1)  //assumes that it is a batch of data
    stdDev/=::sqrt((double)a.d1);
  else
    stdDev/=::sqrt((double)a.N);
  if(!add) for(uint i=0;i<a.N;i++) a.p[i]=stdDev*rnd.gauss();
  else     for(uint i=0;i<a.N;i++) a.p[i]+=stdDev*rnd.gauss();
}

//! a gaussian random vector with Id covariance matrix (sdv = sqrt(dimension))
template<class T>
void rndGauss(MT::Array<T>& a,bool add=false){
  if(!add) for(uint i=0;i<a.N;i++) a.p[i]=rnd.gauss();
  else     for(uint i=0;i<a.N;i++) a.p[i]+=rnd.gauss();
}

//! returns an array with \c dim Gaussian noise elements
template<class T>
MT::Array<T>& rndGauss(double stdDev,uint dim){
  static MT::Array<T> z;
  stdDev/=::sqrt(dim);
  z.resize(dim);
  rndGauss(z,stdDev);
  return z;
}

/*!\brief from a vector of numbers, calculates the softmax distribution
    soft(i) = exp(beta*a(i)), and returns a random sample from
    this distribution (an index in {0,..,a.N-1}) */
template<class T>
uint softMax(const MT::Array<T>& a,arr& soft,double beta){
  double norm=0.,r;
  uint i; int sel=-1;
  soft.resize(a);
  for(i=0;i<a.N;i++){
    soft(i)=exp(beta*a(i));
    norm+=soft(i);
  }
  r=rnd.uni();
  for(i=0;i<a.N;i++){
    soft(i)/=norm;
    r-=soft(i);
    if(sel==-1 && r<0.) sel=i;
  }
  return sel;
}


//===========================================================================
//
//!@name certain initializations
//

/*!\brief return a `dim'-dimensional grid with `steps' intervals
  filling the range [lo,hi] in each dimension */
template<class T>
void grid(MT::Array<T>& a,uint dim,T lo,T hi,uint steps){
  uint i,j;
  if(dim==1){
    a.resize(steps+1,1);
    for(i=0;i<a.d0;i++) a(i,0)=lo+(hi-lo)*i/steps;
    return;
  }
  if(dim==2){
    a.resize(steps+1,steps+1,2);
    for(i=0;i<a.d0;i++) for(j=0;j<a.d1;j++){
      a(i,j,0)=lo+(hi-lo)*j/steps;
      a(i,j,1)=lo+(hi-lo)*i/steps;
    }
    a.reshape(a.d0*a.d1,2);
    return;
  }
  HALT("not implemented yet");
}


//===========================================================================
//
//!@name gnuplot fun
//

//! calls gnuplot to display the (n,2) or (n,3) array (n=number of points of line or surface)
void gnuplot(const arr& X);

//! write 2 arrays in parallel columns in a file
void write(const arr& X,const arr& Y,const char* name);

//! write 3 arrays in parallel columns in a file
void write(const arr& X,const arr& Y,const arr& Z,const char* name);


//===========================================================================
//
//!@name simple image formats
//

/*! save data as ppm or pgm. Images are (height,width,[0,2,3,4])-dim
  byte arrays, where the 3rd dimension determines whether it's a grey
  (0), grey-alpha (2), RGB (3), or RGBA (4) image */
void write_ppm(const byteA &img,const char *file_name,bool swap_rows=false);

/*! read data from an ppm or pgm file */
void read_ppm(byteA &img,const char *file_name,bool swap_rows=false);

//! add an alpha channel to an image array
void add_alpha_channel(byteA &img,byte alpha);

//! make grey scale image
void make_grey(byteA &img);

//===========================================================================
//
//!@name generic data files of double arrays
//

template<class T>
void writeArr(const MT::Array<T>& x,const char* tag,ofstream& os){
  os <<tag <<' ';
  x.writeDim(os);
  os <<'\n';
  os.write((char*)x.p,sizeof(double)*x.N);
  os <<'\n';
}


template<class T>
bool readArr(MT::Array<T>& x,const char *tag,ifstream& is){
  char c;
  if(tag){
    String::readStopSymbols=" \n\r";
    String read_tag;
    is >>read_tag;
    if(!is.good() || read_tag.N()==0) return false;
    CHECK(read_tag==tag,"read `"<<read_tag<<"' instead of `"<<tag<<"' in arr file");
  };
  x.readDim(is);
  c=is.get(); CHECK(c=='\n',"couldn't read newline after ascii tag :-(");
  is.read((char*)x.p,x.sizeT*x.N);
  c=is.get(); CHECK(c=='\n',"couldn't read newline after array buffer :-(");
  return true;
}

inline void scanArrFile(const char* name){
  ifstream is(name,std::ios::binary);
  CHECK(is.good(),"couldn't open file " <<name);
  arr x;
  String tag;
  String::readStopSymbols=" \n\r";
  for(;;){
    is >>tag;
    if(!is.good() || tag.N()==0) return;
    cout <<tag;
    readArr(x,0,is);
    cout <<" =\n" <<x <<endl;
    if(!is.good()) return;
  }
}


//===========================================================================
//
//!@name permutation class
//

namespace MT{
//! A permutation realized as 1-dimensional array of indices.
class Permutation:public Array<uint>{
private:
  Array<uint> storage;

public:
  //! ...
  Permutation(){ init(0); };

public://!@name assignment
  //! sorted permutation of length \c n
  void init(uint n){ resize(n); for(uint i=0;i<N;i++) elem(i)=i; }
  //! random permutation of length \c n
  void random(uint n){ init(n); random(); }
  //! reverse sorted permutation of lenth \c N
  void reverse(uint n){ resize(n); for(uint i=0;i<N;i++) elem(N-1-i)=i; }

public://!@name operations
  //! \c elem(i) = elem(i) + offset modulo permutation range \c N
  void push(int offset=1){ for(uint i=0;i<N;i++) elem(i)=((int)elem(i)+offset)%N; }
  //! store the current permutation
  void store(){ storage.operator=(*this); }
  //! restore the previously stored permutation
  void restore(){ (Array<uint>&)(*this)=storage; }
  //! permute indices of \c i and \c j
  void permute(uint i,uint j){ uint x=elem(i); elem(i)=elem(j); elem(j)=x; }
  //! randomize permutation
  void random(){
    int j,r;
    for(j=N-1;j>=1;j--){
      r=rnd(j+1);
      permute(r,j);
    }
  }
  //! permutes the entries of a according to the current this
  template<class T>
  void permute(Array<T>& a){
    CHECK(N<=a.N,"array smaller than permutation ("<<a.N<<"<"<<N<<")");
    Array<T> b=a;
    for(uint i=0;i<N;i++) a(i)=b(p[i]);
  }
  //! permutes the entries of a according to the inverse of the current this
  template<class T>
  void invpermute(Array<T>& a){
    CHECK(N<=a.N,"array smaller than permutation ("<<a.N<<"<"<<N<<")");
    Array<T> b=a;
    for(uint i=0;i<N;i++) a(p[i])=b(i);
  }
};
}


//===========================================================================
//
// implementations
//

template<class T> char MT::Array<T>::memMoveInit=-1;
template<class T> int MT::Array<T>::sizeT=-1;

#ifdef MT_IMPLEMENTATION
#  include"array.cpp"
#endif


#endif
