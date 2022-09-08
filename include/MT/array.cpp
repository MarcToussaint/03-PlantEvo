/*  Copyright (C) 2000, 2006  Marc Toussaint (mtoussai@inf.ed.ac.uk)
    under the terms of the GNU LGPL (http://www.gnu.org/copyleft/lesser.html)
    see the `std.h' file for a full copyright statement  */

#include "array.h"

//#define MT_CHECK_INVERSE 1e-5

namespace MT{
bool useLapack=true;
}
int ARRAYOLDREAD=0;

static double MT_SIGN_SVD(double a,double b){ return b>0 ? fabs(a) : -fabs(a); }
#define MT_max_SVD(a,b) ( (a)>(b) ? (a) : (b) )

uint own_svd_implementation(const MT::Array<double>& A,
			    MT::Array<double>& U,
			    MT::Array<double>& w,
			    MT::Array<double>& V,
			    bool sort){
  //MT::Array<double*> Apointers,Upointers,Vpointers;
  unsigned m = A.d0; /* rows */
  unsigned n = A.d1; /* cols */
  U.resize(m,n);
  V.resize(n,n);
  w.resize(n);
  double **a = A.getCarray(); //Pointers(Apointers); /* input matrix */
  double **u = U.getCarray(); //Pointers(Upointers); /* left vectors */
  double **v = V.getCarray(); //Pointers(Vpointers); /* right vectors */

  int flag;
  unsigned i,its,j,jj,k,l,nm(0),r;
  double anorm,c,f,g,h,s,scale,x,y,z,t;

  MT::Array<double> rv1(n);

  /* copy A to U */
  for(i=0;i<m;i++) for(j=0;j<n;j++) u[i][j] = a[i][j];

  /* householder reduction to bidiagonal form */
  g = scale = anorm = 0.0;

  for(i=0;i<n;i++){
    l = i + 1;
    rv1(i) = scale * g;
    g = s = scale = 0.0;

    if(i<m){
      for(k=i;k<m;k++) scale += fabs(u[k][i]);

      if(scale!=0.0){
	for(k=i;k<m;k++){
	  u[k][i] /= scale;
	  s += u[k][i] * u[k][i];
	}

	f = u[i][i];
	g = -MT_SIGN_SVD(sqrt(s),f);
	h = f * g - s;
	u[i][i] = f - g;

	for(j=l;j<n;j++){
	  s = 0.0;
	  for(k=i;k<m;k++) s += u[k][i] * u[k][j];

	  f = s / h;
	  for(k=i;k<m;k++) u[k][j] += f * u[k][i];
	}

	for(k=i;k<m;k++) u[k][i] *= scale;
      }
    }

    w(i) = scale * g;
    g = s = scale = 0.0;

    if(i<m && i!=n-1){
      for(k=l;k<n;k++)scale += fabs(u[i][k]);

      if(scale!=0.0){
	for(k=l;k<n;k++){
	  u[i][k] /= scale;
	  s += u[i][k] * u[i][k];
	}

	f = u[i][l];
	g = -MT_SIGN_SVD(sqrt(s),f);
	h = f * g - s;
	u[i][l] = f - g;

	for(k=l;k<n;k++) rv1(k) = u[i][k] / h;

	for(j=l;j<m;j++){
	  s = 0.0;
	  for(k=l;k<n;k++) s += u[j][k] * u[i][k];

	  for(k=l;k<n;k++) u[j][k] += s * rv1(k);
	}

	for(k=l;k<n;k++) u[i][k] *= scale;
      }
    }

    anorm = MT_max_SVD(anorm,fabs(w(i)) + fabs(rv1(i)));
  }

  /* accumulation of right-hand transformations */
  for(l=i=n;i--;l--){
    if(l<n){
      if(g!=0.0){
	  /* double division avoids possible underflow */
	for(j=l;j<n;j++) v[j][i] = (u[i][j] / u[i][l]) / g;

	for(j=l;j<n;j++){
	  s = 0.0;
	  for(k=l;k<n;k++) s += u[i][k] * v[k][j];

	  for(k=l;k<n;k++) v[k][j] += s * v[k][i];
	}
      }

      for(j=l;j<n;j++) v[i][j] = v[j][i] = 0.0;
    }

    v[i][i] = 1.0;
    g = rv1(i);
  }

  /* accumulation of left-hand transformations */
  for(l=i=(m<n?m:n);i--;l--){
    g = w(i);

    for(j=l;j<n;j++) u[i][j] = 0.0;

    if(g!=0.0){
      g = 1.0 / g;

      for(j=l;j<n;j++){
	s = 0.0;
	for(k=l;k<m;k++) s += u[k][i] * u[k][j];

	/* double division avoids possible underflow */
	f = (s / u[i][i]) * g;

	for(k=i;k<m;k++) u[k][j] += f * u[k][i];
      }

      for(j=i;j<m;j++) u[j][i] *= g;
    }else{
      for(j=i;j<m;j++) u[j][i] = 0.0;
    }

    u[i][i]++;
  }

  /* diagonalization of the bidiagonal form */
  for(k=n;k--;){
    for(its=1;its<=30;its++){
      flag = 1;

      /* test for splitting */
      for(l = k + 1;l--;){
	/* rv1 [0] is always zero, so there is no exit */
	nm = l - 1;

	if(fabs(rv1(l)) + anorm == anorm){
	  flag = 0;
	  break;
	}

	//if(!l) break; //(mt 07-01-16)
	if(fabs(w(nm)) + anorm == anorm) break;
      }

      if(flag){
	/* cancellation of rv1 [l] if l greater than 0 */
	c = 0.0;
	s = 1.0;

	for(i=l;i<=k;i++){
	  f = s * rv1(i);
	  rv1(i) *= c;

	  if(fabs(f) + anorm == anorm) break;

	  g = w(i);
	  h = hypot(f,g);
	  w(i) = h;
	  h = 1.0 / h;
	  c = g * h;
	  s = -f * h;

	  for(j=0;j<m;j++){
	    y = u[j][nm];
	    z = u[j][i];
	    u[j][nm] = y * c + z * s;
	    u[j][i] = z * c - y * s;
	  }
	}
      }

      /* test for convergence */
      z = w(k);

      if(l==k){
	if(z<0.0){
	  w(k) = -z;
	  for(j=0;j<n;j++) v[j][k] = -v[j][k];
	}
	break;
      }

      if(its==50) HALT("svd failed");
      //if(its==30) throw k;

      /* shift from bottom 2 by 2 minor */
      x = w(l);
      nm = k - 1;
      y = w(nm);
      g = rv1(nm);
      h = rv1(k);
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = hypot(f,1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + MT_SIGN_SVD(g,f))) - h)) / x;

      /* next qr transformation */
      c = s = 1.0;

      for(j=l;j<k;j++){
	i = j + 1;
	g = rv1(i);
	y = w(i);
	h = s * g;
	g *= c;
	z = hypot(f,h);
	rv1(j) = z;
	c = f / z;
	s = h / z;
	f = x * c + g * s;
	g = g * c - x * s;
	h = y * s;
	y *= c;

	for(jj=0;jj<n;jj++){
	  x = v[jj][j];
	  z = v[jj][i];
	  v[jj][j] = x * c + z * s;
	  v[jj][i] = z * c - x * s;
	}

	z = hypot(f,h);
	w(j) = z;

	/* rotation can be arbitrary if z is zero */
	if(z!=0.0){
	  z = 1.0 / z;
	  c = f * z;
	  s = h * z;
	}

	f = c * g + s * y;
	x = c * y - s * g;

	for(jj=0;jj<m;jj++){
	  y = u[jj][j];
	  z = u[jj][i];
	  u[jj][j] = y * c + z * s;
	  u[jj][i] = z * c - y * s;
	}
      }

      rv1(l) = 0.0;
      rv1(k) = f;
      w(k) = x;
    }
  }

  //sorting:
  if(sort){
    unsigned i,j,k;
    double   p;

    for(i=0;i<n-1;i++){
      p = w(k=i);

      for(j=i+1;j<n;j++) if(w(j)>=p) p = w(k=j);

      if(k!=i){
	w(k) = w(i);
	w(i) = p;

	for(j=0;j<n;j++){
	  p       = v[j][i];
	  v[j][i] = v[j][k];
	  v[j][k] = p;
	}

	for(j=0;j<m;j++){
	  p       = u[j][i];
	  u[j][i] = u[j][k];
	  u[j][k] = p;
	}
      }
    }
  }

  //rank analysis

  for(r=0;r<n && w(r)>0.;r++);

  t = r < n ? fabs(w(n-1)) : 0.0;
  r = 0;
  s = 0.0;
  while(r<n && w(r)>t && w(r)+s>s) s += w(r++);

  return r;
}

uint inverse_SVD(arr& invA,const arr& A){
  unsigned i,j,k,m=A.d0,n=A.d1,r;
  arr U,V,w;
  invA.resize(n,m);
  if(m==0 || n==0) return 0;
  if(m==n && m==1){ invA(0,0)=1./A(0,0); return 0; }
  if(m==n && m==2){ inverse2d(invA,A); return 0; }

  r=svd(A,U,w,V,true);

  for(i=0;i<r;i++) w(i) = 1./w(i);
  for(   ;i<n;i++) w(i) = 0.;

  for(i=0;i<n;i++) for(j=0;j<m;j++){
    double* vi = &V(i,0);
    double* uj = &U(j,0);
    double  t  = 0.;
    for(k=0;k<n;k++) t += vi[k] * w(k) * uj[k];
    invA(i,j) = t;
  }

#ifdef MT_CHECK_INVERSE
  arr D,_D; D.setId(m);
  uint me;
  _D=A*invA;
  double err=maxDiff(_D,D,&me);
  CHECK(err<MT_CHECK_INVERSE ,"inverting failed, error="<<err <<" " <<_D.elem(me) <<"!=" <<D.elem(me) <<"\nA="<<A <<"\ninvA=" <<invA <<"\nA*invA=" <<_D);
#endif

  return r;
}

void pseudoInverse(arr& invA,const arr& A,const arr& W,double eps){
  arr invW,tA,E;
  inverse_SVD(invW,W);
  transpose(tA,A);
  if(!eps){
    invA = invW * tA * inverse(A*invW*tA);
  }else{
    E.setDiag(eps,A.d0);
    invA = invW * tA * inverse(E+A*invW*tA);
  }
}

double _determinant(double **A,uint n){
  if(n==1) return A[0][0];
  if(n==2) return A[0][0]*A[1][1]-A[0][1]*A[1][0];
  uint i,j;
  double d=0;
  double **B=new double*[n-1];
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if(j<i) B[j]=&A[j][1];
      if(j>i) B[j-1]=&A[j][1];
    }
    d+=((i&1)?-1.:1.) * A[i][0] * _determinant(B,n-1);
  }
  delete[] B;
  return d;
}

double determinant(const arr& A){
  CHECK(A.nd==2 && A.d0==A.d1,"determinants require a squared 2D matrix");
  //MT::Array<double*> B;
  A.getCarray(); //Pointers(B);
  return _determinant(A.pp,A.d0);
}

double cofactor(const arr& A,uint i,uint j){
  CHECK(A.nd==2 && A.d0==A.d1,"determinants require a squared 2D matrix");
  arr B=A;
  B.delRow(i);
  B.delColumns(j,1);
  return ((i&1)^(j&1)?-1.:1) * determinant(B);
}

void SUS(const arr& p,uint n,uintA& s){
  //following T. Baeck "EA in Theo. and Prac." p120
  s.resize(n);
  double sum=0,ptr=MT::rnd.uni();
  uint i,j=0;
  for(i=0;i<p.N;i++){
    sum+=p(i)*n;
    while(sum>ptr){ s(j)=i; j++; ptr+=1.; }
  }
  //now, 'sum' should = 'n' and 'ptr' has been 'n'-times increased -> 'j=n'
  CHECK(j==n,"error in rnd::SUS(p,n,s) -> p not normalized?");
}

uint SUS(const arr& p){
  double sum=0,ptr=MT::rnd.uni();
  uint i;
  for(i=0;i<p.N;i++){
    sum+=p(i);
    if(sum>ptr) return i;
  }
  HALT("error in rnd::SUS(p) -> p not normalized? " <<p);
  return 0;
}

void gnuplot(const arr& X){
  if(X.nd==2 && X.d1!=2){ //assume array -> splot
    MT::IOraw=true;
    MT::save(X,"z.pltX");
    gnuplot("splot 'z.pltX' matrix with pm3d,'z.pltX' matrix with lines");
    return;
  }
  if(X.nd==2 && X.d1==2){ //assume curve -> plot
    MT::IOraw=true;
    MT::save(X,"z.pltX");
    gnuplot("plot 'z.pltX' us 1:2");
    return;
  }
  if(X.nd==1){ //assume curve -> plot
    MT::IOraw=true;
    arr Y;
    Y.referTo(X);
    Y.resize(Y.N,1);
    MT::save(Y,"z.pltX");
    gnuplot("plot 'z.pltX' us 1");
    return;
  }
}

void write(const arr& X,const arr& Y,const char* name){
  std::ofstream os;
  MT::open(os,name);
  MT::IOraw=true;
  uint i,j;
  if(X.nd==1){
    for(i=0;i<X.N;i++) os <<X(i) <<' ' <<Y(i) <<std::endl;
  }
  if(X.nd==2){
    for(i=0;i<X.d0;i++){
      for(j=0;j<X[i].N;j++) os <<X[i].elem(j) <<' ';
      for(j=0;j<Y[i].N;j++) os <<Y[i].elem(j) <<' ';
      os <<std::endl;
    }
  }
}

void write(const arr& X,const arr& Y,const arr& Z,const char* name){
  std::ofstream os;
  MT::open(os,name);
  MT::IOraw=true;
  uint i,j;
  if(X.nd==1){
    for(i=0;i<X.N;i++) os <<X(i) <<' ' <<Y(i) <<' ' <<Z(i) <<std::endl;
  }
  if(X.nd==2){
    for(i=0;i<X.d0;i++){
      for(j=0;j<X[i].N;j++) os <<X[i].elem(j) <<' ';
      for(j=0;j<Y[i].N;j++) os <<Y[i].elem(j) <<' ';
      for(j=0;j<Y[i].N;j++) os <<Z[i].elem(j) <<' ';
      os <<std::endl;
    }
  }
}

void write_ppm(const byteA &img,const char *file_name,bool swap_rows){
  CHECK(img.nd==2 || (img.nd==3 && img.d2==3),"only rgb or gray images to ppm");
  ofstream os;
  os.open(file_name, std::ios::out | std::ios::binary);
  if(!os.good()) HALT("could not open file `" <<file_name <<"' for output");
  switch(img.d2){
  case 0:  os <<"P5 "<<img.d1 <<' ' <<img.d0 <<" 255\n";  break; //PGM
  case 3:  os <<"P6 "<<img.d1 <<' ' <<img.d0 <<" 255\n";  break; //PPM
  }
  os.write((char*)img.p,img.N);
}

void read_ppm(byteA &img,const char *file_name,bool swap_rows){
  uint mode, width, height, max;
  ifstream is;
  is.open(file_name, std::ios::in | std::ios::binary);
  if(!is.good()) HALT("could not open file `" <<file_name <<"' for input");
  if(is.get()!='P') HALT("NO PPM FILE:" <<file_name);
  is >>mode;
  if(MT::peerNextChar(is)=='#') MT::skipLine(is);
  is >>width >>height >>max; 
  is.get(); //MUST be a white character if everything went ok
  switch(mode){
  case 5:  img.resize(height,width);    break; //PGM
  case 6:  img.resize(height,width,3);  break; //PPM
  }
  if(!swap_rows){
    is.read((char*)img.p,img.N);
  }else{
    for(uint i=img.d0;i--;) is.read((char*)&img(i,0,0),img.d1*img.d2);
  }
}

void add_alpha_channel(byteA &img,byte alpha){
  uint w=img.d1,h=img.d0;
  img.reshape(h*w,3);
  img.insColumns(3,1);
  for(uint i=0;i<img.d0;i++) img(i,3)=alpha;
  img.reshape(h,w,4);
}

void make_grey(byteA &img){
  CHECK(img.nd==3 && (img.d2==3 || img.d1==4),"makeGray requires color image as input");
  byteA tmp;
  tmp.resize(img.d0,img.d1);
  for(uint i=0;i<img.d0;i++) for(uint j=0;j<img.d1;j++){
    tmp(i,j) = ((uint)img(i,j,0) + img(i,j,1) + img(i,j,2))/3; 
  }
  img=tmp;
}

#ifdef MT_EXPRESSIONS
void assign(arr& x){
  CHECK(x.ex,"self-assignment only if it is an expression");
  MT::Ex *e=x.ex;
  x.init();
  x.ex=e;
  assign(x,x);
  delete x.ex;
  x.ex=0;
}

void assign(arr& x,const arr& a){
  if(!a.ex){ x=a; return; }
  MT::Ex &e=*a.ex;
  if(e.op==MT::UNI){
    arr *A=(arr*)e.A;
    if(A->ex) assign(*A);
    if(!e.trans && e.mul==1 && e.add==0 ){ x=*A; return; }
    if(!e.trans && e.mul==1 ){ scalarPlus(x,*A,*((double*)&e.add)); return; }
    if(!e.trans && e.add==0 ){ scalarMultiplication(x,*A,*((double*)&e.mul)); return; }
    if(e.mul==1 && e.add==0 ){ transpose(x,*A); return; }
    HALT("");
  }else{
    arr *A=(arr*)e.A,*B=(arr*)e.B;
    if(A->ex) assign(*A);
    if(B->ex) assign(*B);
    //bool at,bt;
    //double ac,bc,ap,bp;
    switch(e.op){
    case MT::PROD:
      if(!A->ex && !B->ex){ innerProduct(x,*A,*B); return; }
      HALT("prod");
      break;
    case MT::MUL:
      if(!A->ex && !B->ex){ mult(x,*A,*B); return; }
      HALT("mult");
      break;
    case MT::Div:
      if(!A->ex && !B->ex){ div(x,*A,*B); return; }
      HALT("mult");
      break;
    case MT::OUT:
      if(!A->ex && !B->ex){ outerProduct(x,*A,*B); return; }
      HALT("out");
      break;
    case MT::PLUS:
      if(!A->ex && !B->ex){ plus(x,*A,*B); return; }
      //if(A->ex){ ap=A->ex->add; ac=A->ex->mul; at=A->ex->trans; A=(arr*)A->ex->A; }else{ ap=0; ac=1; at=false; }
      //if(B->ex){ bp=B->ex->add; bc=B->ex->mul; bt=B->ex->trans; B=(arr*)B->ex->A; }else{ bp=0; bc=1; bt=false; }
      //if(!at && !bt && !ap && !bp){ plus(x,ac,*A,bc,*B); return; }
      //if(!at && !bt && !B){ scalarPlus(x,*A,bc); return; }
      HALT("plus");
      break;
    case MT::MINUS:
      if(!A->ex && !B->ex){ minus(x,*A,*B); return; }
      //if(A->ex){ ap=A->ex->add; ac=A->ex->mul; at=A->ex->trans; A=(arr*)A->ex->A; }else{ ap=0; ac=1; at=false; }
      //if(B->ex){ bp=B->ex->add; bc=B->ex->mul; bt=B->ex->trans; B=(arr*)B->ex->A; }else{ bp=0; bc=1; bt=false; }
      //if(!at && !bt && !ap && !bp){ plus(x,ac,*A,-bc,*B); return; }
      //if(!at && !bt && !B){ scalarPlus(x,*A,bc); return; }
      HALT("minus");
      break;
    case MT::UNI:
      HALT("shouldn't be here!");
      break;
    }
    HALT("yet undefined expression");
  }
}
#endif



void getIndexTuple(uintA &I,uint i,const uintA &d){
  uint j;
  CHECK(i<product(d),"out of range");
  I.resize(d.N);
  I.setZero();
  for(j=d.N;j--;){
    I.p[j] = i%d.p[j];
    i -= I.p[j];
    i /= d.p[j];
  }
}

#define DEBUG_TENSOR(x) //x

void tensor(arr &X, const arr &A, const uintA &Aid, const arr &B, const uintA &Bid,uint sum){
  CHECK(A.nd==Aid.N && B.nd==Bid.N,"miss-sized tensor references: " <<A.nd <<"!=" <<Aid.N <<" " <<B.nd <<"!=" <<Bid.N);

  uint n=1+MT::MAX(Aid.max(),Bid.max());
  uint i,j,r,s,N,res;
  intA a(n),b(n);
  uintA d(n),dx(n-sum),I,Ia(A.nd),Ib(B.nd);

  DEBUG_TENSOR(cout <<"Aid=" <<Aid <<" Bid=" <<Bid <<endl;)

  // permutation for A
  a=-1;
  for(i=0;i<A.nd;i++) a(Aid(i))=i;
  //j=A.nd;
  //for(i=0;i<n;i++) if(a(i)==-1){ a(i)=j; j++;  }
  DEBUG_TENSOR(cout <<"permutation for A: " <<a <<endl;)

  //permutation for B
  b=-1;
  for(i=0;i<B.nd;i++) b(Bid(i))=i;
  //j=B.nd;
  //for(i=0;i<n;i++) if(b(i)==-1){ b(i)=j; j++; }
  DEBUG_TENSOR(cout <<"permutation for B: " <<b <<endl;)

  //dimensionalities
  for(i=0;i<n;i++){
    if(a(i)!=-1) r=A.D(a(i)); else r=0;
    if(b(i)!=-1) s=B.D(b(i)); else s=0;
    CHECK(!r || !s || r==s,"inconsistent sharing dimensionalities: " <<r <<"!=" <<s);
    d(i)=MT::MAX(r,s);
  }
  DEBUG_TENSOR(cout <<"full dimensionality d=" <<d <<endl;)

  //total elements:
  N=product(d);
  if(!sum){
    res=1;
    //X.resizeTensor(d);
    CHECK(d==X.d,"for security, please set size before");
  }else{
    dx.resize(d.N-sum);
    res=1;
    for(j=0;j<dx.N;j++) dx(j)=d(j);
    for(   ;j<d .N;j++) res*=d(j);
    //X.resizeTensor(dx);
    CHECK(dx==X.d,"for security, please set size before");
  }
  CHECK(N==X.N*res,"");
  DEBUG_TENSOR(cout <<"dx=" <<dx <<" res=" <<res <<endl; )
  
  //here the copying and multiplying takes place...
  X.setZero();
  for(i=0;i<N;i++){
    getIndexTuple(I,i,d);
    for(j=0;j<A.nd;j++) Ia(j)=I(Aid(j));
    for(j=0;j<B.nd;j++) Ib(j)=I(Bid(j));
    //DEBUG_TENSOR(cout <<"i="<<i <<" I="<<I <<" i/res="<<i/res <<" Ia="<<Ia <<" Ib="<<Ib <<endl;)
    if(!sum){
      X.elem(i) = A(Ia) * B(Ib);
    }else{
      X.elem(i/res) += A(Ia) * B(Ib);
    }
  }
}

void sparseProduct(arr& y,arr& A,const arr& x){
  CHECK(x.nd==1 && A.nd==2 && x.d0==A.d1,"not a proper matrix multiplication");
  if(!A.sparse && !x.sparse){
    innerProduct(y,A,x);
    return;
  }
  if(A.sparse && !x.sparse){
    uint i,j,*k,*kstop;
    y.resize(A.d0); y.setZero();
    double *Ap=A.p;
    uintA *elems=A.sparse;
    for(k=elems->p, kstop=elems->pstop; k!=kstop; Ap++){
      i=*k; k++;
      j=*k; k++;
      y.p[i] += (*Ap) * x.p[j];
    }
    return;
  }
  if(A.sparse && x.sparse){
    uint i,j,n,*k,*kstop,*l,*lstop;
    y.clear(); y.nd=1; y.d0=A.d0; y.sparse=new uintA [2]; y.sparse[1].resize(y.d0); y.sparse[1]=(uint)-1;
    double *xp=x.p;
    uintA *elems,*col;
    elems=x.sparse;
    uint *slot;
    for(k=elems->p, kstop=elems->pstop; k!=kstop; xp++){
      j=*k; k++;
      col=A.sparse+(1+j);
      for(l=col->p, lstop=col->pstop;l!=lstop;){
	i =*l; l++;
	n =*l; l++;
	slot=&y.sparse[1](i);
	if(*slot==(uint)-1){
	  *slot=y.N;
	  y.resizeMEM(y.N+1,true); y(y.N-1)=0.;
	  y.sparse[0].append(i);
	  CHECK(y.sparse[0].N==y.N,"");
	}
	i=*slot;
	y(i) += A.elem(n) * (*xp);
      }
    }
    return;
  }
  if(!A.sparse && x.sparse){
    uint i,j,*k,*kstop,d1=A.d1;
    y.resize(A.d0); y.setZero();
    double *xp=x.p;
    uintA *elems;
    elems=x.sparse;
    for(k=elems->p, kstop=elems->pstop; k!=kstop; xp++){
      j=*k; k++;
      for(i=0; i<A.d0; i++){
	y.p[i] += A.p[i*d1+j] * (*xp);
      }
    }
    return;
  }
}
