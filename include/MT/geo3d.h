/*  Copyright (C) 2000, 2006  Marc Toussaint (mtoussai@inf.ed.ac.uk)
    under the terms of the GNU LGPL (http://www.gnu.org/copyleft/lesser.html)
    see the `std.h' file for a full copyright statement  */

#ifndef MT_geo3d_h
#define MT_geo3d_h

#include "std.h"
#ifdef MT_GL
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif
#include "array.h"

extern void glDrawRect(float,float,float,float,float,float,
		       float,float,float,float,float,float);

extern void glDrawText(const char* txt,float x,float y,float z=0.);


//! A collection of 3D geometry classes
namespace geo3d{
  //forward declarations:
  class Vector;
  class Matrix;
  double operator*(const Vector& a,const Vector& b);
  Vector operator^(const Vector& b,const Vector& c);
  Vector operator*(double b,const Vector& c);
}


namespace geo3d{
//! Simple 3D vector, internally stored as a double[3]
class Vector{
public:
  double v[3]; //<! the vector components

  //! scalar product (inner product)
  friend inline double operator*(const Vector& a,const Vector& b){
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  }

  //! cross product (corresponds to antisymmetric exterior product)
  friend inline Vector operator^(const Vector& b,const Vector& c){
    Vector a;
    a(0)=b[1]*c[2]-b[2]*c[1];
    a(1)=b[2]*c[0]-b[0]*c[2];
    a(2)=b[0]*c[1]-b[1]*c[0];
    return a;
  }

  //! sum of two vectors 
  friend inline Vector operator+(const Vector& b,const Vector& c){
    Vector a;
    a(0)=b[0]+c[0];
    a(1)=b[1]+c[1];
    a(2)=b[2]+c[2];
    return a;
  }

  //! difference between two vectors
  friend inline Vector operator-(const Vector& b,const Vector& c){
    Vector a;
    a(0)=b[0]-c[0];
    a(1)=b[1]-c[1];
    a(2)=b[2]-c[2];
    return a;
  }

  //! multiplication with a scalar
  friend inline Vector operator*(double b,const Vector& c){
    Vector a;
    a(0)=b*c[0];
    a(1)=b*c[1];
    a(2)=b*c[2];
    return a;
  }

  //! multiplication with a scalar
  friend inline Vector operator*(const Vector& b,double c){ return c*b; }

  //! division by a scalar
  friend inline Vector operator/(const Vector& b,double c){ return (1./c)*b; }

  //! multiplication with a scalar
  friend inline Vector& operator*=(Vector& a,double c){
    a(0)*=c; a(1)*=c; a(2)*=c;
    return a;
  }
  //! divide by a scalar
  friend inline Vector& operator/=(Vector& a,double c){
    a(0)/=c; a(1)/=c; a(2)/=c;
    return a;
  }
  //! add a vector
  friend inline Vector& operator+=(Vector& a,const Vector& b){
    a(0)+=b[0]; a(1)+=b[1]; a(2)+=b[2];
    return a;
  }
  //! subtract a vector
  friend inline Vector& operator-=(Vector& a,const Vector& b){
    a(0)-=b[0]; a(1)-=b[1]; a(2)-=b[2];
    return a;
  }
  //! return the negative of a vector
  friend inline Vector operator-(const Vector& b){
    Vector a;
    a(0)=-b[0]; a(1)=-b[1]; a(2)=-b[2];
    return a;
  }

public:
  Vector(){}
  Vector(double x,double y,double z){ set(x,y,z); }

  //! copy operator
  Vector& operator=(const Vector& b){
    v[0]=b[0]; v[1]=b[1]; v[2]=b[2]; return *this;
  }

  //! copy operator
  Vector& operator=(const double* b){
    v[0]=b[0]; v[1]=b[1]; v[2]=b[2]; return *this;
  }

  //! set all entries to same value
  Vector& operator=(double b){
    v[0]=v[1]=v[2]=b; return *this;
  }


public://{ access
  //! lhs reference
  double& operator()(int i){ CHECK(i>=0 && i<3,"geo3d::Vector access - out of range"); return v[i]; }
  double& operator[](int i){ CHECK(i>=0 && i<3,"geo3d::Vector access - out of range"); return v[i]; }

#ifdef MT_MSVC
  //! double-pointer access
  operator double*(){ return v; }
#endif

  //! double-pointer access
  operator const double*() const{ return v; }

  //! set the vector
  void set(double x,double y,double z){ v[0]=x; v[1]=y; v[2]=z; }

  //! set the vector
  void set(double* x){ v[0]=x[0]; v[1]=x[1]; v[2]=x[2]; }

  //! set the vector
  void setZero(){ v[0]=v[1]=v[2]=0.; }

  //! a random vector in [-1,1]^3
  void setRandom(){ v[0]=rnd.uni(-1,1); v[1]=rnd.uni(-1,1); v[2]=rnd.uni(-1,1); }

public://{ vector operations

  //! this=this+b
  void add(double x,double y,double z){ v[0]+=x; v[1]+=y; v[2]+=z; }

  //! this=this-b
  void subtract(double x,double y,double z){ v[0]-=x; v[1]-=y; v[2]-=z; }

  //! this=this/length(this)
  void normalize(){ (*this)/=length(); }

  //! this=this*l/length(this)
  void setLength(double l){
    if(!(*this)) MT_MSG("can't change length of null vector");
    (*this)*=l/length();
  }

  //! this=component of this normal to \c b, (unnormalized!)
  void makeNormal(const Vector& b){
    double l=b.length(),s=v[0]*b[0]+v[1]*b[1]+v[2]*b[2];
    s/=l*l;
    v[0]-=s*b[0]; v[1]-=s*b[1]; v[2]-=s*b[2];
  }

  //! this=component of this colinear to \c b, (unnormalized!)
  void makeColinear(const Vector& b){
    // *this = ((*this)*b)/b.length()) * (*this);
    double l=b.length(),s=v[0]*b[0]+v[1]*b[1]+v[2]*b[2];
    s/=l*l;
    v[0]=s*b[0]; v[1]=s*b[1]; v[2]=s*b[2];
  }

public://{ measuring the vector

  //! is zero?
  bool operator!() const{ return (v[0]==0. && v[1]==0. && v[2]==0.); }

  //! is it normalized?
  bool normalized() const{ return fabs(lengthSqr()-1.)<1e-6; }

  //! returns the length of this
  double length() const{ return ::sqrt(lengthSqr()); }

  //! returns the square of length |a|^2
  double lengthSqr() const{ return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]; }

  //! angle in [0..pi] between this and b
  double angle(const Vector& b) const{
    double x=((*this)*b)/(length()*b.length());
    if(x<-1.) x=-1.;
    if(x>1.) x=1.;
    return ::acos(x);
  }

  /*!\brief if \c this and \c b are colinear, it returns the factor c such
      that this=c*b; otherwise it returns zero */
  double colinear(const Vector& b) const{
    double c=v[0]/b[0];
    if(v[1]==c*b[1] && v[2]==c*b[2]) return c;
    return 0.;
  }


public://{ sphere coordinates

  //! the radius in the x/y-plane
  double radius(){ return ::sqrt(v[0]*v[0]+v[1]*v[1]); }
  //! the angle in the x/y-plane in [-pi, pi]
  double phi(){
    double p;
    if(v[0]==0. || ::fabs(v[0])<1e-10) p=MT_PI/2.; else p=::atan(v[1]/v[0]);
    if(v[0]<0.) if(v[1]<0.) p-=MT_PI; else p+=MT_PI;
    return p;
  }
  //! the angle from the x/y-plane
  double theta(){ return ::atan(v[2]/radius())+MT_PI/2.; }


public://{ I/O
  void write(std::ostream& os) const{
    if(!MT::IOraw) os <<'(' <<v[0] <<' ' <<v[1] <<' ' <<v[2] <<')';
    else os <<' ' <<v[0] <<' ' <<v[1] <<' ' <<v[2];
  }
  void read(std::istream& is){
    if(!MT::IOraw) is >>"(" >>v[0] >>v[1] >>v[2] >>")";
    else is >>v[0] >>v[1] >>v[2];
  }
  //}
};
}
stdPipes(geo3d::Vector);


namespace geo3d{
//! A 3D Matrix, internally stored as double[9]
class Matrix{
public:
  double m[9];

  //! multiplication of two matrices
  friend inline Matrix operator*(const Matrix& b,const Matrix& c){
    Matrix a;
    a(0,0)=b(0,0)*c(0,0)+b(0,1)*c(1,0)+b(0,2)*c(2,0);
    a(0,1)=b(0,0)*c(0,1)+b(0,1)*c(1,1)+b(0,2)*c(2,1);
    a(0,2)=b(0,0)*c(0,2)+b(0,1)*c(1,2)+b(0,2)*c(2,2);

    a(1,0)=b(1,0)*c(0,0)+b(1,1)*c(1,0)+b(1,2)*c(2,0);
    a(1,1)=b(1,0)*c(0,1)+b(1,1)*c(1,1)+b(1,2)*c(2,1);
    a(1,2)=b(1,0)*c(0,2)+b(1,1)*c(1,2)+b(1,2)*c(2,2);

    a(2,0)=b(2,0)*c(0,0)+b(2,1)*c(1,0)+b(2,2)*c(2,0);
    a(2,1)=b(2,0)*c(0,1)+b(2,1)*c(1,1)+b(2,2)*c(2,1);
    a(2,2)=b(2,0)*c(0,2)+b(2,1)*c(1,2)+b(2,2)*c(2,2);
    return a;
  }
  //! sum of two matrices
  friend inline Matrix operator+(const Matrix& b,const Matrix& c){
    Matrix a;
    a.m[0]=b.m[0]+c.m[0]; a.m[1]=b.m[1]+c.m[1]; a.m[2]=b.m[2]+c.m[2];
    a.m[3]=b.m[3]+c.m[3]; a.m[4]=b.m[4]+c.m[4]; a.m[5]=b.m[5]+c.m[5];
    a.m[6]=b.m[6]+c.m[6]; a.m[7]=b.m[7]+c.m[7]; a.m[8]=b.m[8]+c.m[8];
    return a;
  }
  //! transformation of a vector
  friend inline Vector operator*(const Matrix& b,const Vector& c){
    Vector a;
    a(0)=b(0,0)*c.v[0]+b(0,1)*c.v[1]+b(0,2)*c.v[2];
    a(1)=b(1,0)*c.v[0]+b(1,1)*c.v[1]+b(1,2)*c.v[2];
    a(2)=b(2,0)*c.v[0]+b(2,1)*c.v[1]+b(2,2)*c.v[2];
    return a;
  }
  //! multiplication with a scalar
  friend inline Matrix& operator*=(Matrix& a,double c){
    a.m[0]*=c; a.m[1]*=c; a.m[2]*=c;
    a.m[3]*=c; a.m[4]*=c; a.m[5]*=c;
    a.m[6]*=c; a.m[7]*=c; a.m[8]*=c;
    return a;
  }
  //! multiplication with scalar
  friend inline Matrix operator*(double b,const Matrix& c){
    Matrix a;
    a=c;
    a*=b;
    return a;
  }
  //! sum of two matrices
  friend inline Matrix& operator+=(Matrix& a,const Matrix& b){
    a.m[0]+=b.m[0]; a.m[1]+=b.m[1]; a.m[2]+=b.m[2];
    a.m[3]+=b.m[3]; a.m[4]+=b.m[4]; a.m[5]+=b.m[5];
    a.m[6]+=b.m[6]; a.m[7]+=b.m[7]; a.m[8]+=b.m[8];
    return a;
  }

  //! const access via two row and column indices
  const double& operator()(int i,int j) const{ return m[i*3+j]; }
  //! LHS access via two row and column indices
  double& operator()(int i,int j){ return m[i*3+j]; }
private:
  double& M(int i,int j){ return m[i*3+j]; }
public:

  //! reset to zero
  void setZero(){
    m[0]=m[4]=m[8]=
    m[1]=m[2]=m[3]=m[5]=m[6]=m[7]=0.;
  }

  //! reset to identity
  void setId(){
    m[0]=m[4]=m[8]=1.;
    m[1]=m[2]=m[3]=m[5]=m[6]=m[7]=0.;
  }

  //! assign the matrix to the transformation from unit frame to given XYZ frame
  void setFrame(Vector& X,Vector& Y,Vector& Z){
    M(0,0)=X(0); M(0,1)=Y(0); M(0,2)=Z(0);
    M(1,0)=X(1); M(1,1)=Y(1); M(1,2)=Z(1);
    M(2,0)=X(2); M(2,1)=Y(2); M(2,2)=Z(2);
  }
  //! assign the matrix to the transformation from the ORTHOGONAL XYZ frame to the unit frame
  void setInvFrame(Vector& X,Vector& Y,Vector& Z){
    M(0,0)=X(0); M(0,1)=X(1); M(0,2)=X(2);
    M(1,0)=Y(0); M(1,1)=Y(1); M(1,2)=Y(2);
    M(2,0)=Z(0); M(2,1)=Z(1); M(2,2)=Z(2);
  }
  //! assign the matrix to a rotation around the X-axis with angle a (in rad units)
  void setXrot(double a){
    M(0,0)=1.; M(0,1)=0.;     M(0,2)=0.;
    M(1,0)=0.; M(1,1)=cos(a); M(1,2)=-sin(a);
    M(2,0)=0.; M(2,1)=sin(a); M(2,2)= cos(a);
  }

  void setOdeMatrix(double* o){
    m[0]=o[0]; m[1]=o[1]; m[2]=o[2];
    m[3]=o[4]; m[4]=o[5]; m[5]=o[6];
    m[6]=o[8]; m[7]=o[9]; m[8]=o[10];
  }

  void setTensorProduct(const Vector& b,const Vector& c){
    m[0]=b[0]*c[0]; m[1]=b[0]*c[1]; m[2]=b[0]*c[2];
    m[3]=b[1]*c[0]; m[4]=b[1]*c[1]; m[5]=b[1]*c[2];
    m[6]=b[2]*c[0]; m[7]=b[2]*c[1]; m[8]=b[2]*c[2];
  }
};
}


namespace geo3d{
/*!\brief A rotation, internally stored as double[4] quaternion, implements compound,
    multiplication with vector, and export as OpenGL/ODE transformation matrix */
class Rotation{
public:
  double q[4];

  //! initializes to identity
  Rotation(){
    setZero();
  }

  //! compound of two rotations (A=B*C)
  friend inline Rotation operator*(const Rotation& b,const Rotation& c){
    Rotation a;
    a.q[0] = b.q[0]*c.q[0] - b.q[1]*c.q[1] - b.q[2]*c.q[2] - b.q[3]*c.q[3];
    a.q[1] = b.q[0]*c.q[1] + b.q[1]*c.q[0] + b.q[2]*c.q[3] - b.q[3]*c.q[2];
    a.q[2] = b.q[0]*c.q[2] + b.q[2]*c.q[0] + b.q[3]*c.q[1] - b.q[1]*c.q[3];
    a.q[3] = b.q[0]*c.q[3] + b.q[3]*c.q[0] + b.q[1]*c.q[2] - b.q[2]*c.q[1];
    return a;
  }

  //! A=B*C^{-1}
  friend inline Rotation operator/(const Rotation& b,const Rotation& c){
    Rotation a;
    a.q[0] =-b.q[0]*c.q[0] - b.q[1]*c.q[1] - b.q[2]*c.q[2] - b.q[3]*c.q[3];
    a.q[1] = b.q[0]*c.q[1] - b.q[1]*c.q[0] + b.q[2]*c.q[3] - b.q[3]*c.q[2];
    a.q[2] = b.q[0]*c.q[2] - b.q[2]*c.q[0] + b.q[3]*c.q[1] - b.q[1]*c.q[3];
    a.q[3] = b.q[0]*c.q[3] - b.q[3]*c.q[0] + b.q[1]*c.q[2] - b.q[2]*c.q[1];
    return a;
  }

  //! transform of a vector by a rotation
  friend inline Vector operator*(const Rotation& b,const Vector& c){
    double m[9];
    b.getMatrix(m);
    Vector a;
    a.v[0]=m[0]*c.v[0]+m[1]*c.v[1]+m[2]*c.v[2];
    a.v[1]=m[3]*c.v[0]+m[4]*c.v[1]+m[5]*c.v[2];
    a.v[2]=m[6]*c.v[0]+m[7]*c.v[1]+m[8]*c.v[2];
    return a;
  }
  
  //! inverse transform of a vector by a rotation
  friend inline Vector operator/(const Rotation& b,const Vector& c){
    double m[9];
    b.getMatrix(m);
    Vector a;
    a.v[0]=m[0]*c.v[0]+m[3]*c.v[1]+m[6]*c.v[2];
    a.v[1]=m[1]*c.v[0]+m[4]*c.v[1]+m[7]*c.v[2];
    a.v[2]=m[2]*c.v[0]+m[5]*c.v[1]+m[8]*c.v[2];
    return a;
  }
  

public: //{ Manipulations
  //! inverts the current rotation
  void invert(){ q[0]=-q[0]; }

  //! use as similarity measure (distance = 1 - |scalarprod|)
  double scalarProduct(const Rotation& a,const Rotation& b) const{
    return a.q[0]*b.q[0]+a.q[1]*b.q[1]+a.q[2]*b.q[2]+a.q[3]*b.q[3];
  }

  //! multiplies the rotation by a factor f (i.e., makes f-times the rotation)
  void multiply(double f){
    if(q[0]==1. || f==1.) return;
    double phi=acos(q[0]);
    phi*=f;
    q[0]=cos(phi);
    f=sin(phi)/sqrt(q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    q[1]*=f; q[2]*=f; q[3]*=f;
  }

  bool normalized(){
    double n=q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    return fabs(n-1.)<1e-6;
  }

  void normalize(){
    double n=q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    n=sqrt(n);
    q[0]/=n; q[1]/=n; q[2]/=n; q[3]/=n;
  }

  /*!\brief roughly, removes all ``components'' of the rotation that are not
      around the given vector v. More precisely, aligns/projects
      the rotation axis (given by q[1],q[2],q[3] of the quaternion)
      with v and re-normalizes afterwards. */
  void alignWith(const Vector& v){
    double s=q[1]*v[0] + q[2]*v[1] + q[3]*v[2];
    if(!s){ setZero(); return; }// are orthogonal
    s/=v*v;
    q[1]=s*v[0]; q[2]=s*v[1]; q[3]=s*v[2];
    normalize();
  }



public: //{ Assignments
  //! set the quad
  void set(double* x){ q[0]=x[0]; q[1]=x[1]; q[2]=x[2]; q[3]=x[3]; }
  //! set the quad
  void set(double q0,double x,double y,double z){ q[0]=q0; q[1]=x; q[2]=y; q[3]=z; }
  //! reset the rotation to identity
  void setZero(){ q[0]=1; q[1]=q[2]=q[3]=0; }
  //! samples the rotation uniformly from the whole SO(3)
  void setRandom(){
    double s,s1,s2,t1,t2;
    s=rnd.uni();
    s1=sqrt(1-s);
    s2=sqrt(s);
    t1=MT_2PI*rnd.uni();
    t2=MT_2PI*rnd.uni();
    q[0]=cos(t2)*s2;
    q[1]=sin(t1)*s1;
    q[2]=cos(t1)*s1;
    q[3]=sin(t2)*s2;
  }

  //! sets this to a smooth interpolation between two rotations
  void setInterpolate(double t,const Rotation& a,const Rotation b){
    double sign=1.;
    if(a.scalarProduct(a,b)<0) sign=-1.;
    q[0]=a.q[0]+t*(sign*b.q[0]-a.q[0]);
    q[1]=a.q[1]+t*(sign*b.q[1]-a.q[1]);
    q[2]=a.q[2]+t*(sign*b.q[2]-a.q[2]);
    q[3]=a.q[3]+t*(sign*b.q[3]-a.q[3]);
    normalize();
  }

  //! assigns the rotation to \c a DEGREES around the vector (x,y,z)
  void setDeg(double degree,double x,double y,double z){ setRad(degree*MT_PI/180.,x,y,z); }
  void setDeg(double degree,const Vector& vec){ setRad(degree*MT_PI/180.,vec[0],vec[1],vec[2]); }
  //! assigns the rotation to \c a RADIANTS (2*PI-units) around the vector (x,y,z)
  void setRad(double angle,double x,double y,double z){
    double l = x*x + y*y + z*z;
    if(l<1e-15){ setZero(); return; }
    angle/=2.;
    l=sin(angle)/sqrt(l);
    q[0]=cos(angle);
    q[1]=x*l;
    q[2]=y*l;
    q[3]=z*l;
  }
  //! ..
  void setRad(double angle,const Vector &axis){ setRad(angle,axis[0],axis[1],axis[2]); }
  //! assigns the rotation to \c a RADIANTS (2*PI-units) around the current axis
  void setRad(double angle){
    double l = q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    if(l<1e-15){ setZero(); return; }
    angle/=2.;
    l=sin(angle)/sqrt(l);
    q[0]=cos(angle);
    q[1]*=l;
    q[2]*=l;
    q[3]*=l;
  }
  //! sets the rotation as given by a quaternion
  void setQuat(double s,double x,double y,double z){
    q[0]=s; q[1]=x; q[2]=y; q[3]=z;
  }
  //! rotation around the given vector with angle (in rad) equal to norm of the vector
  void setVec(Vector w){
    double phi=w.length();
    setRad(phi,w[0],w[1],w[2]);
  }
  //! rotation that will rotate 'from' to 'to' on direct path
  void setDiff(const Vector& from,const Vector& to){
    double phi=acos(from*to/(from.length()*to.length()));
    Vector axis(from^to);
    if(!axis || !phi) return;
    setRad(phi,axis);
  }


public://{ Access
  //! is identical
  bool operator!() const{ return q[0]==1.; }

#ifdef MT_MSVC
  //! double-pointer access
  operator double*(){ return q; }
#endif

  //! double-pointer access
  operator const double*() const{ return q; }

  //! gets rotation angle (in rad [0,2pi])
  double getRad() const{
    if(q[0]>=1. || q[0]<=-1. || (q[1]==0. && q[2]==0. && q[3]==0.)) return 0;
    return 2.*acos(q[0]);
  }

  //! gets rotation angle (in degree [0,360])
  double getDeg() const{
    if(q[0]>=1. || q[0]<=-1. || (q[1]==0. && q[2]==0. && q[3]==0.)) return 0;
    return 360./MT_PI*acos(q[0]);
  }

  //! gets rotation angle (in degree [0,360]) and vector
  void getDeg(double& degree,Vector& vec) const{
    if(q[0]>=1. || q[0]<=-1. || (q[1]==0. && q[2]==0. && q[3]==0.)){ degree=0.; vec.set(0.,0.,1.); return; }
    degree=acos(q[0]);
    double s=sin(degree);
    degree*=360./MT_PI;
    vec(0)=q[1]/s; vec(1)=q[2]/s; vec(2)=q[3]/s;
  }

  //! gets rotation angle (in rad [0,2pi]) and vector
  void getRad(double& angle,Vector& vec) const{
    if(q[0]>=1. || q[0]<=-1. || (q[1]==0. && q[2]==0. && q[3]==0.)){ angle=0.; vec.set(0.,0.,1.); return; }
    angle=acos(q[0]);
    double s=sin(angle);
    angle*=2;
    vec(0)=q[1]/s; vec(1)=q[2]/s; vec(2)=q[3]/s;
    CHECK(angle>=0. && angle<=MT_2PI,"");
  }

  //! gets the axis rotation vector with length equal to the rotation angle in rad
  Vector& getVec(Vector& v) const{
    if(q[0]>=1. || q[0]<=-1. || (q[1]==0. && q[2]==0. && q[3]==0.)){ v.setZero(); return v; }
    double phi=acos(q[0]);
    double s=2.*phi/sin(phi);
    v(0)=s*q[1]; v(1)=s*q[2]; v(2)=s*q[3];
    return v;
  }

  void getX(Vector& Rx) const{ Rx = (*this)*Vector(1,0,0); }
  void getY(Vector& Ry) const{ Ry = (*this)*Vector(0,1,0); }
  void getZ(Vector& Rz) const{ Rz = (*this)*Vector(0,0,1); }

  //! exports the rotation to a double[9] matrix, row-by-row
  double* getMatrix(double* m) const{
    double qq1 = 2*q[1]*q[1];
    double qq2 = 2*q[2]*q[2];
    double qq3 = 2*q[3]*q[3];
    m[0]=1-qq2-qq3;               m[1]=2*(q[1]*q[2]-q[0]*q[3]); m[2]=2*(q[1]*q[3] + q[0]*q[2]);
    m[3]=2*(q[1]*q[2]+q[0]*q[3]); m[4]=1-qq1-qq3;               m[5]=2*(q[2]*q[3] - q[0]*q[1]);
    m[6]=2*(q[1]*q[3]-q[0]*q[2]); m[7]=2*(q[2]*q[3]+q[0]*q[1]); m[8]=1-qq1-qq2;
    return m;
  }

  void setMatrix(double* m){
    q[0]=sqrt(1.-(3.-(m[0]+m[4]+m[8]))/4.);
    q[3] = (m[3]-m[1])/(4.*q[0]);
    q[2] = (m[2]-m[6])/(4.*q[0]);
    q[1] = (m[7]-m[5])/(4.*q[0]);
    normalize();
    //CHECK(normalized(),"failed :-(");
  }

  //! exports the rotation to an ODE format matrix of type double[12]
  double* getOdeMatrix(double* m) const{
    double qq1 = 2*q[1]*q[1];
    double qq2 = 2*q[2]*q[2];
    double qq3 = 2*q[3]*q[3];
    m[0]=1-qq2-qq3;               m[1]=2*(q[1]*q[2]-q[0]*q[3]); m[2] =2*(q[1]*q[3] + q[0]*q[2]);
    m[4]=2*(q[1]*q[2]+q[0]*q[3]); m[5]=1-qq1-qq3;               m[6] =2*(q[2]*q[3] - q[0]*q[1]);
    m[8]=2*(q[1]*q[3]-q[0]*q[2]); m[9]=2*(q[2]*q[3]+q[0]*q[1]); m[10]=1-qq1-qq2;
    m[3]=m[7]=m[11]=0.;
    return m;
  }
#ifdef MT_GL
  //! exports the rotation to an OpenGL format matrix of type GLfloat[16]
  GLfloat* getOpenGLMatrix(GLfloat* m) const{
    GLfloat qq1 = 2*q[1]*q[1];
    GLfloat qq2 = 2*q[2]*q[2];
    GLfloat qq3 = 2*q[3]*q[3];
    m[0]= 1 - qq2 - qq3;
    m[4] = 2*(q[1]*q[2] - q[0]*q[3]);
    m[8] = 2*(q[1]*q[3] + q[0]*q[2]);
    m[1] = 2*(q[1]*q[2] + q[0]*q[3]);
    m[5] = 1 - qq1 - qq3;
    m[9] = 2*(q[2]*q[3] - q[0]*q[1]);
    m[2] = 2*(q[1]*q[3] - q[0]*q[2]);
    m[6] = 2*(q[2]*q[3] + q[0]*q[1]);
    m[10] = 1 - qq1 - qq2;
    m[3]=m[7]=m[11]=m[12]=m[13]=m[14]=0.;
    m[15]=1.;
    return m;
  }
#endif

public://{ I/O
  void writeNice(std::ostream& os) const{ Vector v; os <<"Rotation: " <<getDeg() <<" around " <<getVec(v) <<"\n"; }
  void write(std::ostream& os) const{
    if(!MT::IOraw) os <<'(' <<q[0] <<' ' <<q[1] <<' ' <<q[2] <<' ' <<q[3] <<')';
    else os <<' ' <<q[0] <<' ' <<q[1] <<' ' <<q[2] <<' ' <<q[3];
  }
  void read(std::istream& is){ is >>"(" >>q[0] >>q[1] >>q[2]  >>q[3] >>")"; }
  //}
};
}
stdPipes(geo3d::Rotation);


namespace geo3d{
/*! 3D inertial, ortho-normal frame. A frame can equally be used to represent
    an absolute coordinate frame (`a frame relative to world coordinates') as well as
    a transformation between frames (`a frame relative to another frame') */
class Frame{
public:
  Vector p;  //!< position
  Rotation r;//!< orientation (rotation w.r.t. world frame)
  Vector v;  //!< linear velocity
  Vector w;  //!< angular velocity
  Vector a;  //!< acceleration
  Vector b;  //!< angular acceleration
  double s;  //!< scale

public://{ basics
  //! constructor
  Frame(){ reset(); }

  //! copy operator
  Frame& operator=(const Frame& f){
    p=f.p; v=f.v; r=f.r; w=f.w; a=f.a; b=f.b; s=f.s; return *this; }

  //! copy operator
  Frame& set(const Frame& f){
    p=f.p; v=f.v; r=f.r; w=f.w; a=f.a; b=f.b; s=f.s; return *this; }

  //! initialize by reading from the string
  Frame& set(const char* txt){ set(MT::String(txt)()); return *this; }

  //! initialize by reading from the stream
  Frame& set(istream &is){ reset(); read(is); return *this; }

  //! resets the position to origin, rotation to identity, velocities to zero, scale to unit
  void reset(){
    p.setZero(); v.setZero(); r.setZero(); w.setZero(); a.setZero(); b.setZero(); s=1.; }

  //! randomize the frame
  void random(){
    p.setRandom(); v.setRandom(); r.setRandom(); w.setRandom(); a.setRandom(); b.setRandom(); s=rnd.uni(); }

  /*!\brief moves the frame according to the current velocities \c v and \c w
      and the time span \c time (if time is given in seconds, v has
      dimension units/sec, and w has dimension rad/sec) */
  void step(double time=1.){
    Rotation W;
    W.setVec(w);
    W.multiply(time);
    p+=v;
    r=W*r;
  }

public://{ turtle commands -- transformations RELATIVE to the current frame
  //! multiply the current scale by f
  void scale(double f){
    s*=f;
  }
  //! move the turtle by the vector (x,z,y) WITH RESPECT TO the current orientation/scale
  void translate(double x,double y,double z){
    Vector X(x,y,z);
    X=r*(s*X); //in global coords
    p+=X;
    v+=w^X;
  }
  //! add a velocity to the turtle's inertial frame
  void addVelocity(double x,double y,double z){
    Vector X(x,y,z);
    v+=r*(s*X);
  }
  //! rotate the turtle orientation by an angle (given in DEGREE) around the vector (x,y,z) (given relative to the current orientation)
  void rotateDeg(double degree,double x,double y,double z){
    Rotation R;
    R.setDeg(degree,x,y,z);
    r=r*R;
  }
  //! rotate the turtle orientation by an angle (given in radiants) around the vector (x,y,z) (given relative to the current orientation)
  void rotateRad(double rad,double x,double y,double z){
    Rotation R;
    R.setRad(rad,x,y,z);
    r=r*R;
  }
  //! rotate the turtle orientation as given by a quaternion
  void rotateQuat(double s,double x,double y,double z){
    Rotation R;
    R.q[0]=s; R.q[1]=x; R.q[2]=y; R.q[3]=z;
    r=r*R;
  }
  //! add an angular velocity to the turtle inertial frame
  void addAngularVelocityDeg(double degree,double x,double y,double z){
    Vector W(x,y,z); W.normalize();
    W*=degree*MT_PI/180.;
    w+=r*W;
  }
  //! add an angular velocity to the turtle inertial frame
  void addAngularVelocityRad(double rad,double x,double y,double z){
    Vector W(x,y,z); W.normalize();
    W*=rad;
    w+=r*W;
  }
  //! add an angular velocity to the turtle inertial frame
  void addAngularVelocityRad(double wx,double wy,double wz){
    Vector W(wx,wy,wz);
    w+=r*W;
  }
  /*!\brief transform the turtle into the frame f,
      which is interpreted RELATIVE to the current frame
      (new = f * old) */
  void transform(const Frame& f){
    Vector P(r*(s*f.p)); //relative offset in global coords
    Vector V(r*(s*f.v)); //relative vel in global coords
    Vector W(r*f.w);     //relative ang vel in global coords
    p += P;
    v += w^P;
    v += V;
    a += b^P;
    a += w^((w^P) + 2.*V);
    a += r*(s*f.a);
    b += w^W;
    b += r*f.b;
    w += W;
    r = r*f.r;
    s*=f.s;
  }
  //! inverse transform (new = f^{-1} * old) or (old = f * new) 
  void invtransform(const Frame& f){
    s/=f.s;
    r=r/f.r;
    w-=r*f.w;
    v-=r*(s*f.v);
    Vector P(r*(s*f.p)); //frame offset in global coords
    v-=w^P;
    p-=P;
  }
  //! new = old * f
  void makeRelative(const Frame& f){
    Frame newf;
    newf.set(f);
    newf.transform(*this);
    *this=newf;
  }
  //! new = old * f^{-1}
  void makeRelativeToInv(const Frame& f){
    Frame newf;
    newf.setInverse(f);
    newf.transform(*this);
    *this=newf;
  }
  //! this = f^{-1}
  void setInverse(const Frame& f){ reset(); invtransform(f); }
  //!  to = new * from
  void setDifference(const Frame& from,const Frame& to){
    s=to.s/from.s;
    r=Rotation()/from.r *to.r;
    w=from.r/(to.w-from.w);
    v=(1./from.s) * (from.r/(to.v-from.v));
    v-=(1./from.s) * (from.r/(from.w^(to.p-from.p)));
    p=(1./from.s) * (from.r/(to.p-from.p));
  }

public://{ export the frame
#ifdef MT_GL
  //! get the current position/orientation/scale in an OpenGL format matrix (of type GLfloat[16])
  GLfloat* getOpenGLMatrix(GLfloat *m){
    double M[9]; r.getMatrix(M);
    m[0]=s*M[0]; m[4]=s*M[1]; m[8] =s*M[2]; m[12]=p[0];
    m[1]=s*M[3]; m[5]=s*M[4]; m[9] =s*M[5]; m[13]=p[1];
    m[2]=s*M[6]; m[6]=s*M[7]; m[10]=s*M[8]; m[14]=p[2];
    m[3]=0.;     m[7]=0.;     m[11]=0.;     m[15]=1.;
    return m;
  }

  //! get inverse OpenGL matrix for this frame (of type GLfloat[16]) */
  GLfloat* getOpenGLInverseMatrix(GLfloat *m){
    double M[9]; r.getMatrix(M);
    Vector pinv; pinv=r/p;
    m[0]=s*M[0]; m[4]=s*M[3]; m[8] =s*M[6]; m[12]=-pinv[0];
    m[1]=s*M[1]; m[5]=s*M[4]; m[9] =s*M[7]; m[13]=-pinv[1];
    m[2]=s*M[2]; m[6]=s*M[5]; m[10]=s*M[8]; m[14]=-pinv[2];
    m[3]=0.;     m[7]=0.;     m[11]=0.;     m[15]=1.;
    return m;
  }

  void glLoadMatrix(){
    GLfloat GLmatrix[16];
    getOpenGLMatrix(GLmatrix);
    glLoadMatrixf(GLmatrix);
  }
#endif

public://{ I/O
  //! operator<<
  void write(std::ostream& os) const{
    if(!!p) os <<" t" <<p;
    if(!!v) os <<" v" <<v;
    if(!!w) os <<" w" <<w;
    if(!!r) os <<" q" <<r;
    if(s!=1.) os <<" s(" <<s <<") ";
    if(!!p || !!v || !!r || !!w || s!=1.) os <<' ';
  }
  //! operator>>
  void read(std::istream& is){
    reset();
    char c;
    double x[4];
    for(;;){
      is >>c;
      if(is.fail()) return; //EOF I guess
      //if(c==';') break;
      //if(c==',') is >>c;
      switch(c){
      case 't': is>>"(">>x[0]>>x[1]>>x[2]>>")";	      translate(x[0],x[1],x[2]); break;
      case 'q': is>>"(">>x[0]>>x[1]>>x[2]>>x[3]>>")"; rotateQuat(x[0],x[1],x[2],x[3]); break;
      case 'r': is>>"(">>x[0]>>x[1]>>x[2]>>x[3]>>")"; rotateRad(x[0],x[1],x[2],x[3]); break;
      case 'd': is>>"(">>x[0]>>x[1]>>x[2]>>x[3]>>")"; rotateDeg(x[0],x[1],x[2],x[3]); break;
      case 'v': is>>"(">>x[0]>>x[1]>>x[2]>>")";	      addVelocity(x[0],x[1],x[2]); break;
      case 'w': is>>"(">>x[0]>>x[1]>>x[2]>>")";       addAngularVelocityRad(x[0],x[1],x[2]); break;
      case 's': is>>"(">>x[0]>>")";                   scale(x[0]); break;
      default: is.putback(c); return; //all other symbols finish the reading
      }
      if(is.fail()) HALT("error reading `"<<c<<"' parameters in frame");
    }
    if(is.fail()) HALT ("could not read Frame struct");
  }
  //}
};
}
stdPipes(geo3d::Frame);



namespace geo3d{
/*!\brief A turtle is basically an orthogonal frame (defined here as position, orientation, and scale)
  which can be moved around by turtle commands (which are always relative to the current position/orientation/scale)
  and which implements a push/pop stack of frames */
class Turtle:public Frame{
public:
  MT::Array<Frame> stack;

public:
  Turtle(){ reset(); stack.memMove=true; }


public://{ turtle commands
  //! remember the current turtle position/orientation/scale to be recalled later
  void push(){ stack.append(*this); }
  //! reinstall the last turtle position/orientation/scale and remove it from memory
  void pop(){ repop(); stack.remove(stack.N-1); }
  //! reinstall the last turtle position/orientation/scale, but keep in memory to recall again
  void repop(){
    CHECK(stack.N>=1,"no frame to pop");
    ((Frame&)(*this))=stack(stack.N-1);
  }
};
}


namespace geo3d{
//! camera in OpenGL convention
class Camera{
public:
  Frame *f;
  Vector foc;
  bool ownFrame;
  Vector offset;
  float phi,theta;

  float heightAbs;
  float heightAngle;
  float whRatio;
  float zNear,zFar;

  /*!\brief constructor; specify a frame if the camera is to be attached
      to an existing frame. Otherwise the camera creates its own
      frame */
  Camera(Frame *frame=0){
    if(!frame){ f=new Frame; ownFrame=true; }else{ f=frame; ownFrame=false; }
    reset();
  }
  ~Camera(){
    if(ownFrame) delete f;
  }

  void reset(){
    f->reset();
    offset.setZero();
    theta=phi=0.;
    heightAngle=90.;
    heightAbs=10.;
    whRatio=1.;
    zNear=1.;
    zFar=1000.;
  }

public://{ set camera's perspective:
  //! the height angle (in degrees) of the camera perspective; set it 0 for orthogonal projection
  void setHeightAngle(float a){ heightAngle=a; }
  //! the absolute height of the camera perspective (automatically also sets heightAngle=0)
  void setHeightAbs(float h){ heightAngle=0; heightAbs=h; }
  //! the z-range (depth range) visible for the camera
  void setZRange(float znear,float zfar){ zNear=znear; zFar=zfar; }
  //! set the width/height ratio of your viewport to see a non-distorted picture
  void setWHRatio(float ratio){ whRatio=ratio; }
  //! the frame's position
  void setPosition(float x,float y,float z){ f->p.set(x,y,z); }
  //! the camera's additional offset within the frame
  void setOffset(float x,float y,float z){ offset.set(x,y,z); }
  //! rotate the frame to focus the absolute coordinate origin (0,0,0)
  void focusOrigin(){ foc.setZero(); focus(); }
  //! rotate the frame to focus the point (x,y,z)
  void focus(float x,float y,float z){ foc.set(x,y,z); focus(); }
  //! rotate the frame to focus the point given by the vector
  void focus(const Vector& v){ foc=v; focus(); }
  //! rotate the frame to focus (again) the previously given focus
  void focus(){ Vector v(foc-f->p); watchDirection(v[0],v[1],v[2]); } //f->Z=f->p; f->Z-=foc; f->Z.normalize(); upright(); }
  //! rotate the frame to watch in the direction vector D
  void watchDirection(float x,float y,float z){
#if 1
    Vector D(x,y,z),V(0,0,-1);
    f->r.setZero();
    if(D[0]==0. && D[1]==0.){
      if(D[2]>0) f->r.setDeg(180,1,0,0);
      return;
    }
    D.normalize();
    V=f->r*V;
    Rotation r;
    r.setDiff(V,D);
    f->r=r*f->r;
    upright();
#else
    f->Z.set(-D[0],-D[1],-D[2]); f->Z.normalize();
    upright();
#endif
  }
  //! rotate the frame to set it upright (i.e. camera's y aligned with 's z)
  void upright(){
#if 1
    //construct desired X:
    Vector V(0,0,-1),X(1,0,0),dX,up;
    X=f->r*X; //true X
    V=f->r*V;
    if(fabs(V[2])<1.) up.set(0,0,1); else up.set(0,1,0);
    dX=up^V; //desired X
    if(dX*X<=0) dX=-dX;
    Rotation r;
    r.setDiff(X,dX);
    f->r=r*f->r;
#else
    if(f->Z[2]<1.) f->Y.set(0,0,1); else f->Y.set(0,1,0);
    f->X=f->Y^f->Z; f->X.normalize();
    f->Y=f->Z^f->X; f->Y.normalize();
#endif
  }

  //}

#ifdef MT_GL
  /*! sets OpenGL's GL_PROJECTION matrix accordingly -- should be
      called in an opengl draw routine */
  void glSetProjectionMatrix(){
    if(heightAngle==0)
      glOrtho(-whRatio*heightAbs/2,whRatio*heightAbs/2,
	      -heightAbs/2,heightAbs/2,zNear,zFar);
    else
      gluPerspective(heightAngle,whRatio,zNear,zFar);
    glRotatef(theta,1,0,0);
    glRotatef(phi,0,1,0);
    glTranslatef(-offset[0],-offset[1],-offset[2]);
    GLfloat m[16];
    glMultMatrixf(f->getOpenGLInverseMatrix(m));
  }

  //! convert from gluPerspective's non-linear [0,1] depth to the true [zNear,zFar] depth
  void glConvertToTrueDepth(double &d){
    d = zNear + (zFar-zNear)*d/(zFar/zNear*(1.-d)+1.);
  }

  //! convert from gluPerspective's non-linear [0,1] depth to the linear [0,1] depth
  void glConvertToLinearDepth(double &d){
    d = d/(zFar/zNear*(1.-d)+1.);
  }
#endif

  //! a GL draw routine that draws a little red cone to indicate the camera
  void glDraw(){
    HALT("not implemented yet!");
  }
  static void staticDraw(void* classP){ ((Camera*)classP)->glDraw(); }
};
}


namespace geo3d{
#define MT_FIG_header "#FIG 3.2\nLandscape\nCenter\nInches\nLetter\n100.00\nSingle\n-2\n1200 2"
//! a 3D path
class Path{
public:
  MT::Array<float> a;
  bool clip;
  Path(){ clip=false; }

  //! number of points
  uint N(){ return a.d0; }
  //! delete all points
  void clear(){ a.resize(0,3); }
  //! append a point
  void append(float x,float y,float z){
    uint n=N();
    a.resizeCopy(n+1,3);
    a(n,0)=x; a(n,1)=y; a(n,2)=z;
  }
  //! saves the path in an xfig file
  void saveFIG(const char* name){
    std::ofstream os(name);
    CHECK(os.good() && !os.bad(),"could not open file '" <<name <<"' for output");
    uint i;
    os <<MT_FIG_header <<std::endl;
    os <<"0 32 #000000\n";
    if(!clip){
      os <<"2 1 0 1 32 0 50 -1 -1 0.000 0 0 -1 0 0 " <<N() <<"\n";
      for(i=0;i<N();i++)
        os <<"\t" <<(int)floor(5000.*a(i,0)) <<" " <<-(int)floor(5000.*a(i,1)) <<"\n";
    }else{
      for(i=0;i+1<N();i+=2){
        os <<"2 1 0 1 32 0 50 -1 -1 0.000 0 0 -1 0 0 2\n";
        os <<"\t" <<(int)floor(5000.*a(i,0)) <<" " <<-(int)floor(5000.*a(i,1)) <<"\n";
        os <<"\t" <<(int)floor(5000.*a(i+1,0)) <<" " <<-(int)floor(5000.*a(i+1,1)) <<"\n";
      }
    }
    os.close();
  }
  static void staticDraw(void* classP){ ((Path*)classP)->glDraw(); }
  void glDraw(){
#ifdef MT_GL
    uint i;
    double s=.02;
    glColor3f(0,0,0);
    for(i=0;i<N();i++){
      glDrawRect(-s+a(i,0),-s+a(i,1),.01+a(i,2),
		  s+a(i,0),-s+a(i,1),.01+a(i,2),
		  s+a(i,0), s+a(i,1),.01+a(i,2),
		 -s+a(i,0), s+a(i,1),.01+a(i,2));
      if(i){
	glBegin(GL_LINES);
	glVertex3f(a(i-1,0),a(i-1,1),a(i-1,2));
	glVertex3f(a(i  ,0),a(i  ,1),a(i  ,2));
	glEnd();
      }
      glColor3f(0,0,0.7);
      glDrawText(STRING(i),a(i,0)+s,a(i,1)+s,a(i,2)+.05);
    }
#endif
  }
};
}
//stdPipes(Path);


namespace geo3d{
/*!\brief a 3D triangularized shape. Composed of a set of vertices, a set of
  triangles (each stored as point index triplets), a set of strips,
  and additional properties like normals and colors for vertices or
  triangles */
class Trimesh{
public:
  arr V;                //!< vertices
  uintA T;              //!< triangles (faces)
  arr N;                //!< vertex normals
  MT::Array<Color> C;   //!< vertex colors
  MT::Array<uintA> strips; //!< triangle strips (each with a 1D stripe index set)

  intA G;               //!< group for each vertex
  MT::Array<Frame*> GF; //!< frame for each group (GF.N is number of groups)
  MT::Array<uintA>  GT; //!< triangles for each group (GT.N=GF.N+1, last entry contains mixed group triangles)
  bool colored;

  Trimesh(){ colored=false; }

private:
  void firstPass(FILE* file);
  void secondPass(FILE* file); 

public:
  uint &Tni(uint v,uint i){ static uint dummy; return dummy; } //normal index
  uint &Tti(uint v,uint i){ NIY; return T(v,i); } //texture index

public:
  /*!\ initialises the ascii-obj file "filename" */
  void readOBJ(char* filename);

  void clear();
  /*!\brief calculate the normals of all triangles (Tn) and the average
    normals of the vertices (Vn); average normals are averaged over
    all adjacent triangles that are in the triangle list or member of
    a strip */
  void calcNormals();

  /*!\brief add triangles according to the given grid; grid has to be a 2D
    Array, the elements of which are indices referring to vertices in
    the vertex list (V) */
  void gridToTriangles(const MT::Array<uint> &grid);

  /*!\brief add strips according to the given grid (sliced in strips along
    the x-axis (the first index)); grid has to be a 2D Array, the
    elements of which are indices referring to vertices in the vertex
    list (V) */
  void gridToStrips(const MT::Array<uint>& grid);

  /*!\brief add strips according to the given grid (sliced in strips along
    the x-axis (the first index)); it is assumed that the vertices in
    the list V linearly correspond to points in the XxY grid */
  void gridToStrips(uint X,uint Y);

  /*!\brief add triangles according to the given grid; grid has to be a 2D
    Array, the elements of which are indices referring to vertices in
    the vertex list (V) */
  void gridToTriangles(uint X,uint Y);

  /*!\brief delete all void triangles (with vertex indices (0,0,0)) and void
    vertices (not used for triangles or strips) */
  void deleteUnused();

  //!\brief internal draw routine for OpenGL
  void glDraw();

  void makeRel(){
    uint i;
    int g;
    Vector *v;
    for(i=0;i<V.d0;i++) if((g=G(i))!=-1){
      v = (Vector*)&V(i,0);
      *v = GF(g)->r/((*v) - GF(g)->p);
      v = (Vector*)&N(i,0);
      *v = GF(g)->r/(*v);
    }
  }
  void collectTriGroups(){
    uint i;
    int g;
    GT.resize(GF.N+1);
    for(i=0;i<T.d0;i++){
      g=G(T(i,0));
      if(g!=-1 && g==G(T(i,1)) && g==G(T(i,2))){
	GT(g).append(i);
      }else{
	GT(GF.N).append(i);
      }
    }
  }
};
}

#ifdef MT_IMPLEMENTATION
#  include"geo3d.cpp"
#endif

#endif
