/*  Copyright (C) 2000, 2006  Marc Toussaint.
    email: mtoussai@inf.ed.ac.uk

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    A copy of the GNU Lesser General Public License can usually be found
    at http://www.gnu.org/copyleft/lesser.html; if not, write to the Free
    Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
    02110-1301 USA */

/*! \file std.h
    \brief My base file: contains a bunch of standard includes, macros, operators,
    etc. that are used by most other headers */
#ifndef MT_std_h
#define MT_std_h


//--------- I need a system flag!
#if !defined MT_MSVC && !defined MT_Cygwin && !defined MT_Linux && !defined MT_MinGW
#  error please define a system flag (one of MT_MSVC, MT_Linux, MT_Cygwin, or MT_MinGW)! (for MSVC also define /Zm500)
#endif

 #include <string.h>

//--------- standard includes:
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
//#include <typeinfo>
#if defined MT_MSVC // || defined WIN32 || defined _WIN32 || defined _MSC_VER
//#  include<limits>
#  include<time.h>
#  include<sys/timeb.h> 
#  define MT_TIMEB
#  ifdef MT_QT
#    undef  NOUNICODE
#    define NOUNICODE
#  endif
#  pragma warning(disable: 4305 4244 4250 4355 4786)
#endif
#if defined MT_MinGW
#  include<unistd.h>
#  include<sys/time.h>
#  include<sys/timeb.h>
#  define MT_TIMEB
#endif
#if defined MT_Linux || defined MT_Cygwin
//#  include<unistd.h>
#  include<sys/time.h>
#  include<sys/times.h>
#  include<sys/resource.h>
#endif

#if defined MT_QT || defined MT_FREEGLUT
#  ifdef MT_MSVC
#    include<windows.h>
#    undef min //I hate it that windows defines these macros!
#    undef max
#  endif
#  define MT_GL
#  ifdef MT_QT
#    undef MT_FREEGLUT
#  endif
#endif
#ifdef MT_doxy
#  define DOX !\brief
#endif

#ifndef MT_QT
#  define Q_OBJECT
#endif

//--------- declare my namespace for the first time:
//! Marc Toussaint namespace
namespace MT{
  inline int debugBreakpoint(){
    char c; c=0; //insert a break point here, if you want...
    return 0; 
  }
}


//--------- error handling:
#ifndef HALT
#  define MT_MSG(msg) { std::cerr <<"MSG  @ " <<__FILE__ <<':' <<__LINE__\
          <<":\n --- "<<msg<<" --- "<<std::endl; MT::debugBreakpoint(); }
#  define HALT(msg) {   std::cerr <<"HALT @ " <<__FILE__ <<':' <<__LINE__\
          <<":\n --- "<<msg<<" --- "<<std::endl; MT::debugBreakpoint(); abort(); }
#  define NIY HALT("not implemented yet")
#endif


//--------- check macros:
#ifndef MT_NOCHECK
#  define CHECK(cond,msg) if(!(cond)) HALT("CHECK failed: "<<msg);
#  define IFCHECK(x) x
#else
#  define CHECK(cond,msg)
#  define IFCHECK(x)
#endif


//--------- basic defs:
#define MT_PI 3.14159265358979323846
//#define MT_2PI (2.*3.14159265358979323846)
#define MT_LN2 0.69314718055994528622676398299518041312694549560546875
#define MT_2PI 6.283195307179587
#define MT_LnSqrt2Pi -0.9189385332046727417803296 
#define MT_SQRT2 1.414213562373095049
#define MT_SQRTPI 1.772453850905516027
typedef unsigned char byte;      //!< byte
typedef unsigned short int word; //!< 2 bytes
typedef unsigned int uint;       //!< unsigned integer
typedef enum _FwdBwd { FORWARD, BACKWARD } FwdBwd;


//--------- macros for looping
#ifndef forall
#define forall(i,a)            for(a.First(i)              ;a.CheckNext(i);    a.Next(i))
#define forall_rev(i,a)        for(a.Last(i)               ;a.CheckPrev(i);    a.Prev(i))
#define forall_save(i,j,a)     for(a.First(i),j=i,a.Next(j);a.CheckNext(i);i=j,a.Next(j))
#define forall_rev_save(i,j,a) for(a.Last(i) ,j=i,a.Prev(j);a.CheckPrev(i);i=j,a.Prev(j))
#endif


//--------- macros to define the standard << and >> operatos for most my classes:
#define stdInPipe(type)\
  inline std::istream& operator>>(std::istream& is,type& x){ x.read(is);return is; }
#define stdEqPipe(type)\
  inline type& operator<<(type& x,const char* str){ std::istringstream ss(str); ss >>x; return x; }
#define stdOutPipe(type)\
  inline std::ostream& operator<<(std::ostream& os,const type& x){ x.write(os); return os; }
#define stdPipes(type)\
  inline std::istream& operator>>(std::istream& is,type& x){ x.read(is);return is; }\
  inline std::ostream& operator<<(std::ostream& os,const type& x){ x.write(os); return os; }\
//! assignes the object x to the description given by str
#define fwdStdPipes(type)\
  inline std::istream& operator>>(std::istream& is,type& x);\
  inline std::ostream& operator<<(std::ostream& os,const type& x);


//--------- macros for piping doubles EXACTLY (without rounding errors) in hex coding:
#define OUTHEX(y) "0x" <<std::hex <<*((unsigned long*)&y) <<std::dec
#define INHEX(y) std::hex >>*((unsigned long*)&y) >>std::dec


//--------- tools to scan (parse) from a stream"
namespace MT{
  //! skip all kinds of white characters in the stream
  inline char skip(std::istream&,char* skipchars=" \n\r\t");

  //! tell you about the next char (after skip()) but puts it back in the stream
  inline char peerNextChar(std::istream& is){ char c; skip(is); is.get(c); is.putback(c); return c; }
}
//! a global operator to scan (parse) strings from a stream
inline std::istream& operator>>(std::istream& is,const char* str){
  CHECK(is.good(),"bad enter");
  uint i,n=strlen(str);
  char *buf=new char [n+1]; buf[n]=0;
  MT::skip(is," \n\r\t");
  is.read(buf,n);
  if(!is.good() || strcmp(str,buf)){
    for(i=n;i--;) is.putback(buf[i]);
    is.setstate(std::ios::failbit);
    std::cout <<"WARNING: scanning of constant string `" <<str
	      <<"' failed! (read instead: `" <<buf <<"')" <<std::endl;
  }
  delete[] buf;
  return is;
}
//! a global operator to scan (parse) strings from a stream
/*inline std::istream& operator>>(std::istream& is,const char c){
  char d;
  MT::skip(is," \n\r\t");
  is.get(d);
  if(!is.good() || c!=d){
    is.putback(d);
    is.setstate(std::ios::failbit);
    std::cout <<"WARNING: scanning of constant char `" <<c
	      <<"' failed! (read instead: `" <<d <<"')" <<std::endl;
  }
  return is;
}*/
//! the same global operator for non-const string
inline std::istream& operator>>(std::istream& is,char* str){
  return is >>(const char*)str;
}


//--------- standard conversion base type to a class:
#ifndef MT_doxy
//! the basetype template `casts' an elemental type into a class
template<class T>
class basetype{
public:
  T x; //!< the single datum
  basetype():x(){};
  basetype(const T _x):x(_x){};

  T& data(){ return x; }
  operator T() const{ return x; } //!< implicit and const access
  void read(std::istream& is){ is >>x; } //!< implements standard operator>>
  void write(std::ostream& os) const{ os <<x; } //!< implements standard operator<<
};
template<class T> stdInPipe(basetype<T>);
template<class T> stdOutPipe(basetype<T>);

namespace MT{
  typedef basetype<byte> Byte;
  typedef basetype<int> Int;
  typedef basetype<float> Float;
  typedef basetype<double> Double;
}
#endif



//===========================================================================
//
// a very basic memory buffer
//

namespace MT{
//! minimal memory buffer using memove internally; a simplification of the Mem class
template<class T>
class Mem{
public:
  T *p;    //!< pointer to memory
  uint N;  //!< \# elements
  uint M;  //!< actual buffer size (in terms of # elements)
  uint S;  //!< element size (sizeof(T))
  Mem(){ p=0; S=sizeof(T); N=0; M=0; }
  Mem(uint n){ p=0; S=sizeof(T); N=0; M=0; resize(n); }
  ~Mem(){ if(M) delete[] p; }
  void append(const T& x){ resize(N+1); operator()(N-1)=x; }
  void remove(uint i,uint n=1){
    if(N>i+n) memmove(p+i,p+i+n,S*(N-i-n));
    resizeCopy(N-n);
  }
  void free(){ if(M) delete[] p; p=0; M=N=0; }
  void zero(){ memset(p,0,S*M); }
  T &operator()(uint i) const{ CHECK(i<N,"Mem range error ("<<i<<"<"<<N<<")"); return p[i]; }
  void resizeCopy(uint n){ resize(n,true); }
  void resize(uint n,bool copy=true){
    if(N==n) return;
    T *pold=p;
#if 0 //rigid allocation: M always equals exactly N
    if(M!=n){
      if(n) p=new T [n];
      if(!p) HALT("MT::Mem failed memory allocation of "<<n*S<<"bytes");
      if(copy) memmove(p,pold,S*(N<n?N:n));
      if(M) delete[] pold;
      N=M=n;
    }
#else //flexible allocation (more than needed in case of multiple resizes)
    if(n>0 && M==0){ //first time
      p=new T [n];
      if(!p) HALT("MT::Mem failed memory allocation of "<<n*S<<"bytes");
      M=n;
    }else if(n>M || 10+2*n<M/2){
      uint oversize=10+2*n;
      p=new T [oversize];
      if(!p) HALT("MT::Mem failed memory allocation of "<<n*S<<"bytes");
      if(copy) memmove(p,pold,S*(N<n?N:n));
      if(M) delete[] pold;
      M=oversize;
    }
    N=n;
#endif
  }
  Mem& operator=(const Mem& m){
    CHECK(S==m.S,"");
    resize(m.N);
    memmove(p,m.p,S*N);
    return *this;
  }
};
}

typedef MT::Mem<double> doubleM;
typedef MT::Mem<float>  floatM;
typedef MT::Mem<uint>   uintM;
typedef MT::Mem<int>    intM;
typedef MT::Mem<char>   charM;
typedef MT::Mem<byte>   byteM;
typedef MT::Mem<bool>   boolM;


//===========================================================================
//
// standard little methods in my namespace
//

namespace MT{
  extern int argc;
  extern char** argv;
  extern std::ifstream cfgFile;
  extern bool cfgOpenFlag,cfgLock;
  extern Mem<void*> parameters;
  //! stream modifier for some classes (Mem in particular)
  extern bool IOraw;
  //! no logfile: default=true, becomes false when MT::init is called
  extern bool noLog;

  //! access to the log-file
  std::ofstream& log(const char *name="MT.log");

  //! open an output-file with name '\c name'
  void open(std::ofstream& fs,const char *name,char* errmsg="");

  //! open an input-file with name '\c name'
  void open(std::ifstream& fs,const char *name,char* errmsg="");

  /*!\brief a standard method to save an object into a file. The same as
    std::ofstream file; MT::open(file,filename); file <<x;
    file.close(); */
  template<class T>
  void save(const T& x,const char *filename){
    std::ofstream file;
    MT::open(file,filename);
    file <<x;
    file.close();
  }

  /*!\brief a standard method to load object from a file. The same as
    std::ifstream file; MT::open(file,filename); file >>x;
    file.close(); */
  template<class T>
  void load(T& x,const char *filename){
    std::ifstream file;
    MT::open(file,filename);
    file >>x;
    file.close();
  }

  //! returns true if the (0-terminated) string s contains c 
  bool contains(const char* s,char c);

  //! skips the chars (typically white characters) when parsing from the istream
  char skip(std::istream& is,char* skipchars);

  //! skips a newline character (same as skip(is,"\n");)
  void skipLine(std::istream& is);

  //! returns the i-th of str 
  byte bit(byte *str,uint i);

  //! flips the i-th bit of b
  void flip(byte& b,uint i);

  //! filps the i-th bit of b
  void flip(int& b,uint i);

  /*!\brief the distance between x and y w.r.t.\ a circular topology
      (e.g. modMetric(1,8,10)=3) */
  double modMetric(double x,double y,double mod);

  //! the sign (+/-1) of x (+1 for zero)
  double sign(double x);

  //! returns 0 for x<0, 1 for x>1, x for 0<x<1
  double linsig(double x);

  void constrain(double& x,double a,double b);

  //! the angle of the vector (dx,dy) in [-pi,pi]
  double phi(double dx,double dy);
  double dphi(double x,double y,double dx,double dy);

  /*!\brief save division, checks for division by zero; force=true will return
    zero if y=0 */
  double DIV(double x,double y,bool force=false);

  inline double sigmoid11(double x){ return x/(1.+::fabs(x)); }
  inline double sigmoid(double x){ return 1./(1.+exp(-x)); }
  inline double dsigmoid(double x){ double y=sigmoid(x); return y*(1.-y); }

  template<class T> T MIN(T a,T b){ return a<b?a:b; }
  template<class T> T MAX(T a,T b){ return a>b?a:b; }

  //! approximate exp (sets up a static value table)
  double approxExp(double x);

  //! ordinary Log, but cutting off for small values
  double Log(double x);

  //! integer log2
  uint Log2(uint n);

  //! square of a double
  double sqr(double x);

  //!
  double sinc(double x);

  //!
  double cosc(double x);

  /*!\brief real time since start of the process in floating-point seconds
    (probably in micro second resolution) -- Windows checked! */
  double realTime();

  /*!\brief user CPU time of this process in floating-point seconds (pure
    processor time) -- Windows checked! */
  double cpuTime();

  /*!\brief system CPU time of this process in floating-point seconds (the
    time spend for file input/output, x-server stuff, etc.)
    -- not implemented for Windows! */
  double sysTime();

  /*!\brief total CPU time of this process in floating-point seconds (same
    as cpuTime + sysTime) -- not implemented for Windows! */
  double totalTime();
  
  //! the absolute real time and date as string
  char *date();

  //! wait real time
  void wait(double sec);

  //! wait for an ENTER at the console
  void wait();

  //! the integral shared memory size -- not implemented for Windows!
  long mem();

  //! push reset on the stop watch (user CPU time)
  void resetTimer();

  //! read the stop watch (user CPU time)
  double getTimer(bool reset=true,double relTime=-1.);

  //! memorize the command line arguments and open a log file
  void init(int _argc, char *_argv[]);

  /*!\brief Open a (possibly new) config file with name '\c name'.<br> If
    \c name is not specified, it searches for a command line-option
    '-cfg' and, if not found, it assumes \c name=MT.cfg */
  void openConfigFile(char* name=0);

  //! returns false if command-line option \c -tag was not found
  bool checkOption(const char *tag);

  //! returns false if command-line option \c -tag was not found
  bool getOption(const char *tag,char* &option);

  /*!\brief Search for a command line option \c -tag and, if found, pipe the
     next command line option into \c value by the
     \c operator>>(istream&,type&). Returns false on failure. */
  template<class T>
  bool getFromCmdLine(T& x,const char *tag){
    char* opt;
    if(!getOption(tag,opt)) return false;
    std::istringstream s(opt);
    s >>x;
    if(s.fail()) HALT("error when reading parameter from command line: " <<tag);
    return true;
  }

  /*!\brief Search the first occurence of a sequence '\c tag:'
  in the config file (opened automatically) and, if found, pipes
  it in \c value. Returns false if parameter is not found. */
  template<class T>
  bool getFromCfgFile(T& x,const char *tag){
    if(!cfgOpenFlag) openConfigFile();
    CHECK(!cfgLock,"cfg file is locked");
    cfgLock=true;
    cfgFile.clear();
    cfgFile.seekg(std::ios::beg);
    if(!cfgFile.good()){ cfgLock=false; return false; }
    
    unsigned n=strlen(tag);
    char *buf=new char [n+2]; memset(buf,0,n+2);
    while(cfgFile.good()){
      memmove(buf,buf+1,n);
      buf[n]=cfgFile.get();
      if(buf[n]==' ' || buf[n]=='\t' || buf[n]==':' || buf[n]=='='){ buf[n]=0; if(!strcmp(tag,buf)) break; buf[n]=':'; }
    };
    delete[] buf;
    
    if(!cfgFile.good()){ cfgLock=false; return false; }
    
    skip(cfgFile," :=\n\r\t");
    cfgFile >>x;
    
    if(cfgFile.fail()) HALT("error when reading parameter " <<tag);
    cfgLock=false;
    return true;
  }
  
  template<class T>
  void getParameterBase(T& x,const char *tag,bool hasDefault,const T* Default){
    if(!cfgOpenFlag) openConfigFile();
    log() <<std::setw(20) <<tag <<" = " <<std::setw(5);
    log().flush();
    
    if(getFromCmdLine(x,tag)){
      log() <<x <<" [" <<typeid(x).name() <<"] (cmd line!)" <<std::endl;
      return;
    }
    
    if(getFromCfgFile(x,tag)){
      log() <<x <<" [" <<typeid(x).name() <<"]" <<std::endl;
      return;
    }
    
    if(hasDefault){
      x=*Default;
      log() <<x <<" [" <<typeid(x).name() <<"] (default!)"<<std::endl;
      return;
    }
    
    HALT("could not initialize parameter `" <<tag
	     <<"': parameter has no default;\n     either use command option `-"
	     <<tag <<" ...' or specify `"
	     <<tag <<"= ...' in the config file `MT.cfg'");
  }

  template<class T>
  void getParameter(T& x,const char *tag,const T& Default){ getParameterBase(x,tag,true,&Default); }
  template<class T>
  void getParameter(T& x,const char *tag){ getParameterBase(x,tag,false,(T*)0); }

#ifndef MT_doxy
#ifdef MT_QT
  struct InitQt{ InitQt(); static char *className(){ return 0; } static void staticMetaObject(){} };
#else
  struct InitQt{};
#endif
#endif
}



//===========================================================================
//
// Parameter class - I use it frequently to read parameters from file or cmd line
// 

namespace MT{
  /*!\brief A parameter that initializes itself from the command line
    (use \c MT::init), parameter file, or a default value (priority in
    this order).  Initialization is done on the fly the _first_ time
    its value is queried (i.e., referenced by the cast operators).*/
  template<class type>
  class Parameter{
  public:
    const char* typeName;
    type value,Default;
    char* tag;
    bool initialized,hasDefault;

  public:
    //!@name constructors

    //! Determines the tag to search for in parameter file/command line
    Parameter(char* _tag){
      typeName=typeid(type).name();
      MT::parameters.append(this);
      initialized=false;
      tag=_tag;
      hasDefault=false;
    };

    /*!\brief specifies also a default value -- parameter does not have to but
      can be specified in the parameter file/command line */
    Parameter(char* _tag,const type& _default){
      typeName=typeid(type).name();
      MT::parameters.append(this);
      initialized=false;
      tag=_tag;
      hasDefault=true;
      Default=_default;
    };

    ~Parameter(){
      for(uint i=0;i<MT::parameters.N;i++){
	if(MT::parameters(i)==this){ MT::parameters.remove(i); return; }
      }
    }

    //!@name value access

    //! standard type conversion: returns a const of the parameter value
    operator type(){ if(!initialized) initialize(); return value; }

    //! ()-operator: returns an lvalue of the parameter value
    type& operator()(){ if(!initialized) initialize(); return value; }


    //!@name manipulation

    //! assigs a value to the parameter -- no further initialization needed
    type& operator=(const type v){ initialized=true; value=v; return value; }

    //! set the tag (replacing the one from the constructor)
    void setTag(char* _tag){ tag=_tag; }

    /*!\brief enforces that the parameter is reinitialized from the parameter
      file/command line, the next time it is referenced -- even if it
      has been initialized before */
    void reInitialize(){ initialized=false; }


    //!@name explicit grabbing


  private:
    void initialize(){
      if(!initialized){
	getParameterBase(value,tag,hasDefault,&Default);
	initialized = true;
      }
    }
  };

}



//===========================================================================
//
// random number generator
//

namespace MT{
  /*!\brief A random number generator. An global instantiation \c
    MT::rnd of a \c Rnd object is created. Use this one object to get
    random numbers.*/
  class Rnd{
  private:
    bool ready;
    int rpoint;     /* Feldindex    */
    long rfield[256];   /* Schieberegisterfeld  */


  public:
    //! ...
    Rnd(){ ready=false; };


  public://!@name initialization
    //! initialize with a specific seed
    unsigned long seed(unsigned long n);

    //! use Parameter<uint>("seed") as seed
    unsigned long seed();

    //! uses the internal clock to generate a seed
    unsigned long clockSeed();

  public://!@name access
    //! a initeger random number uniformly distributed in [0,?]
    inline uint num(){ if(!ready) seed(); return (uint)rnd250() >> 5; }
    //! same as \c num()
    inline uint operator ()(){ return num(); }
    //! a initeger random number uniformly distributed in [0,\c i-1]
    inline uint num(uint limit){
      CHECK(limit,"zero limit in rnd.num()"); return num() % limit; }
    inline uint num(int lo,int hi){ return lo+num(hi-lo+1); }
    //! same as \c num(i)
    inline uint operator ()(uint i){ return num(i); }
    inline uint operator ()(int lo,int hi){ return num(lo,hi); }
    //! a random variable uniformly distributed in [0,1]
    inline double uni(){ return ((double)num(1<<22))/(1<<22); }
    //! a random variable uniformly distributed in [\c low, \c high]
    inline double uni(double low,double high){ return low+uni()*(high-low); }
    //! a gaussian random variable with mean zero
    double gauss();
    /*!\brief a positive integer drawn from a poisson distribution with given
      \c mean; is case \c mean>100, a (positive) gauss number
      \c floor(mean+gauss(sqrt(mean))+.5) is returned */
    uint poisson(double mean);
    //@}


  private:
    inline long  rnd250(){
      rpoint = (rpoint+1) & 255;          // Index erhoehen
      return rfield[rpoint] =  rfield[(rpoint-250) & 255]
	                     ^ rfield[(rpoint-103) & 255];
    }

    void seed250(long seed);
  };

  //! The global Rnd object
  extern Rnd rnd;
}



//===========================================================================
//
// string class
//

namespace MT{
  typedef std::iostream IOStream;
  
  /*!\brief String implements the functionalities of an ostream and an
      istream, but also can be send to an ostream or read from an
      istream. It is based on a simple streambuf derived from the
      MT::Mem class */
  class String:public IOStream{
  private:
    class StringBuf:public std::streambuf{
    public:
      Mem<char> mem;
      const char* ptr() const{ return mem.p; }
      char* ptr(){ return mem.p; }
      virtual int overflow(int C = traits_type::eof()){
        mem(mem.N-1)=C; mem.append(0); resetI();
        return C;
      }
      void resetI(){ setg(ptr(),ptr(),ptr()+mem.N-1); }
      char *getIpos(){ return gptr(); }
      void setIpos(char* p){
        setg(ptr(),p,ptr()+mem.N-1);
      }
    };
    
  public:
    StringBuf buf;
    static char* readSkipSymbols;
    static char* readStopSymbols;
    static bool eatStopSymbol;
    //! standard constructor
    String():IOStream(&buf){ clr(); }
    //! copy constructor
    String(const String& s):IOStream(&buf){ this->operator=(s); }
    //! copy constructor for an ordinary C-string (needs to be 0-terminated)
    String(const char* s):IOStream(&buf){ this->operator=(s); }
    //! returns a reference to this
    IOStream& str(){ return (IOStream&)(*this); }
    /*!\brief returns the true memory buffer (C-string) of this class (which is
      always kept 0-terminated) */
    operator char*(){ return buf.ptr(); }
    //! as above but const
    operator const char*() const{ return buf.ptr(); }
    //! also returns the char-pointer
    char* ptr(){ return buf.ptr(); }
    //! returns the i-th char
    char operator()(uint i){ return buf.mem(i); }
    String& operator()(){ return *this; }
    //! copy operator
    void operator=(const String& s){ buf.mem=s.buf.mem; buf.resetI(); }
    //! copies from the C-string
    void operator=(const char *s){ buf.mem.resize(strlen(s)+1); memmove(buf.ptr(),s,strlen(s)+1); buf.resetI(); }
    //! shorthand for the !strcmp command
    bool operator==(const char *s){ return !strcmp(buf.ptr(),s); }
    //! shorthand for the !strcmp command
    bool operator==(const String& s){ return !strcmp(buf.ptr(),s.buf.ptr()); }
    //! deletes all memory and resets all stream flags
    String& clr(){ buf.mem.resize(1); buf.mem(0)=0; clear(); return *this; }
    /*!\brief when using this String as an istream (to read other variables
      from it), this method resets the reading-pointer to the beginning
      of the string and also clears all flags of the stream */
    String& resetI(){ buf.resetI(); clear(); return *this; }
    char *getIpos(){ return buf.getIpos(); }
    void setIpos(char* p){ buf.setIpos(p); }
    
    //! returns the length of the string (excluding the terminating 0)
    uint N() const{ return buf.mem.N-1; }
    uint length() const{ return buf.mem.N-1; } // renamed (HH)
    /*!\brief appends some other text/variable/whatever to the string (as if
      piped into it) */
    template<class T>
    String operator+(const T& v) const{
      String news(*this); news <<v; return news;
    }
    //! writes the string into some ostream
    void write(std::ostream& os) const{ os <<buf.ptr(); }
    /*!\brief reads the string from some istream until one of the stopSymbols
      is encountered (default: newline symbols) */
    void read(std::istream& is,char *stopSymbols=NULL){
      if(!stopSymbols) stopSymbols=readStopSymbols;
      MT::skip(is,readSkipSymbols);
      buf.mem.resize(10);
      uint n=0;
      char c=is.get();
      while(is.good() && !MT::contains(stopSymbols,c)){
	if(buf.mem.N<=n) buf.mem.resize(n+10);
	buf.mem(n)=c;
	n++;
	c=is.get();
      }
      if(!eatStopSymbol) is.putback(c);
      is.clear();
      buf.mem.resize(n+1);
      buf.mem(n)=0;
      resetI();
    }
  };
}
stdPipes(MT::String);

#define STRING(x) ((MT::String&)(MT::String().str()<<x)).ptr()



//===========================================================================
//
// color class
//

#ifndef MT_doxy
namespace MT{
  //! simple float[3] color class
  class Color{
  public:
    float
      r, //!< red
      g, //!< green
      b; //!< blue
    
    //! ...
    friend inline Color operator+(const Color& c1,const Color& c2){
      return Color(c1.r+c2.r, c1.g+c2.g, c1.b+c2.b); }
    
    //! ...
    friend inline Color operator*(float f,const Color& c2){
      return Color(f*c2.r, f*c2.g, f*c2.b); }
    
  public:
    //! initializes to white
    Color(){ setGray(1.); }
    
    //! initialize with standard color (using setIndex)
    Color(uint i){ setIndex(i); }

    //! initialize with RGB
    Color(float red,float green,float blue){ setRgb(red,green,blue); }
    
    //! copy operator
    Color& operator=(const Color& c){ r=c.r; g=c.g; b=c.b; return *this; }
    
    //! return true iff black
    bool operator!(){ if(r==0. && g==0. && b==0.) return true; return false; }
    
    //! float-pointer access
    operator const float*() const{ return (float*)this; }
    
    //! chooses color from a color table (distributed around the hue-scale)
    void setIndex(unsigned i){
      if(!i) setRgb(0.,0.,0.); else setHsv(((i-1)*63)%360,255,255); }
    
    //! set RGA values
    void setRgb(float red,float green,float blue){ r=red; g=green; b=blue; }
    
    //! set RGA values from bytes in [0,255]
    void setRgbByte(byte red,byte green,byte blue){
      r=red/255.f; g=green/255.f; b=blue/255.f; }
    
    //! set color by hue [0,360], saturation [0,255], and value [0,255]
    void setHsv(int hue,byte sat,byte val){
      float h=hue/60.f,s=sat/255.f,v=val/255.f;
      h=(float)fmod(h,6.f);
      r=g=b=0.;
      if(h<=1.)        { r=v; g=v*h; }
      if(h>1. && h<=2.){ g=v; r=v*(2.f-h); }
      if(h>2. && h<=3.){ g=v; b=v*(h-2.f); }
      if(h>3. && h<=4.){ b=v; g=v*(4.f-h); }
      if(h>4. && h<=5.){ b=v; r=v*(h-4.f); }
      if(h>5. && h<=6.){ r=v; b=v*(6.f-h); }
      r=s*r+(1.f-s)*v;
      g=s*g+(1.f-s)*v;
      b=s*b+(1.f-s)*v;
    }
    
    //! set color by temperature: hot=red, middle=yellow, cold=blue
    void setTemp(float temp){
      Color hot(1.,0.,0.),middle(1.,1.,0.),cold(0.,0.,1.);
      if(temp>1.) temp=1.;
      if(temp<0.) temp=0.;
      if(temp>.5){ temp=2.f*temp-1.f; *this=temp*hot + (1.-temp)*middle; }
      else{ temp=2.f*temp; *this=temp*middle + (1.f-temp)*cold; }
    }
    
    //! set color by temperature: red - yellow - gray(middle) - green - blue
    void setTemp2(float temp){
      Color r(1.,0.,0.),y(1.,1.,0.),zero(.5,.5,.5),g(0.,1.,0.),b(0.,0.,1.);
      if(temp>1.) temp=1.;
      if(temp<-1.) temp=-1.;
      if(temp>.5){  temp=2.*temp-1.; *this=temp*r + (1.-temp)*y; return; }
      if(temp>.0){  temp=2.*temp;    *this=temp*y + (1.-temp)*zero; return; }
      if(temp>-.5){ temp=-2.*temp;   *this=temp*g + (1.-temp)*zero; return; }
      { temp=-2.*temp-1.;*this=temp*b + (1.-temp)*g; return; }
    }
    
    //! set gray value [0,1]
    void setGray(float gray){ if(gray<0) gray=0.; if(gray>1) gray=1.; r=g=b=gray; }
    
    //! get RGB values as bytes [0,255]
    void getRgb(byte& R,byte& G,byte& B) const{
      R=(byte)(255.*r); G=(byte)(255.*g); B=(byte)(255.*b); }
    
    //! get the gray value (average of RGB)
    float getGray() const{ return (r+g+b)/3.; }
    
    //! mix with white
    void whiten(float f){
      if(f>1.) f=1.; else if(f<0.) f=0.;
      r+=f*(1.-r); g+=f*(1.-g); b+=f*(1.-b);
    }
    
    //! mix with black
    void blacken(float f){
      if(f>1.) f=1.; else if(f<0.) f=0.;
      r-=f*r; g-=f*g; b-=f*b;
    }
    
    //! plain color mixing
    void mix(Color& A,Color& B,float f=.5){
      if(f>1.) f=1.; else if(f<0.) f=0.;
      r=f*A.r+(1.-f)*B.r;
      g=f*A.g+(1.-f)*B.g;
      b=f*A.b+(1.-f)*B.b;
    }
    
    //! additive color mixing
    void mixAdd(Color& A,Color& B,float f=.5){
      if(f>1.) f=1.; else if(f<0.) f=0.;
      r=1.-f*(1.-A.r)+(1.-f)*(1.-B.r);
      g=1.-f*(1.-A.g)+(1.-f)*(1.-B.g);
      b=1.-f*(1.-A.b)+(1.-f)*(1.-B.b);
    }
    
    //! subtractive color mixing
    void mixSub(Color& A,Color& B,float f=.5){
      if(f>1.) f=1.; else if(f<0.) f=0.;
      r=1.-::pow(1.f-A.r,f)*::pow(1.f-B.r,1.f-f);
      g=1.-::pow(1.f-A.g,f)*::pow(1.f-B.g,1.f-f);
      b=1.-::pow(1.f-A.b,f)*::pow(1.f-B.b,1.f-f);
    }
    
    //! take smaller of the two values
    void min(Color& A,Color& B){
      r=A.r<B.r?A.r:B.r;
      g=A.g<B.g?A.g:B.g;
      b=A.b<B.b?A.b:B.b;
    }

    //! prototype for operator <<
    void write(std::ostream& os) const{ os <<"(" <<r <<":" <<g <<":" <<b <<")"; }
    
    //! prototype for operator >>
    void read(std::istream& is){ is >>"(" >>r >>":" >>g >>":" >>b >>")"; }
  };
}
stdPipes(MT::Color);
#endif


//===========================================================================
//
// Bag - a generic container: list of any-types
//

#define BAG_DEBUG(x)

namespace MT{
  /*!\brief a simple any-type container: basically a list of items,
    where an item it a type-identifying tag and a void pointer */
  class Bag{
  public:
    class Item{
    public:
      Item *next; 
      void *value;
      int type;
      Item(){ next=0; value=0; type=-1; }
      virtual ~Item(){}
      virtual void write(std::ostream &os) = 0;
      virtual void free() = 0;
      virtual void set(void *x) = 0;
      //virtual get(void *x) = 0;
    };
    template<class T>
    class TypedItem:public Item{
    public:
      virtual void write(std::ostream &os){
	if(!value) return;
	if(!ItemTypes(type).n){
	  os <<*((T*)value);
	}else{
	  T *t=(T*)value;
	  uint n=ItemTypes(type).n;
	  for(uint i=0;i<n;i++,t++) os <<' ' <<*t;
	}
      }
      virtual void free(){
	if(!value) return;
        if(!ItemTypes(type).n){
	  delete ((T*)value);
	}else{
	  delete[] ((T*)value);
	}
	value=0;
      }
      virtual void set(void *_x){
	if(!_x){ value=0; return; }
	uint n=ItemTypes(type).n;
	T *t=(T*)value;
	T *x=(T*)_x;
	if(value) free();
	if(!n){
	  value = new T;
	  t=(T*)value;
	  *t=*x;
	}else{
	  value = new T[n];
	  t=(T*)value;
	  for(uint i=0;i<n;i++) t[i]=x[i];
	}
      }
      //virtual get(void *x){ *((T*)x) = *((T*)value); }
    };
    struct ItemType{ char *tag; const char *typeidname; uint n; };
    
    Item *first; //!< pointer to first item in the bag
    static Mem<ItemType> ItemTypes; //!< table of types shared by _all_ Bags

    Bag(){ first=0; }

    //! add an arbitrary-type item to the bag
    template<class T>
    void add(const char* tag,T *x,uint n=0){
      Item *a=new TypedItem<T>;
      a->type=getType(tag,typeid(T).name(),n);
      a->set(x);
      a->next=first;
      first=a;
    }

    //! delete an item from the bag (previously retrieved via get)
    void del(Item *a){
      Item *b;
      while(a){
	b=a->next;
	a->free();
	delete a;
	a=b;
      }
    }

    //! write the bag to a stream (prototype for operator<<)
    void write(std::ostream &os) const{
      uint i=0;
      Item *a;
      for(a=first;a;a=a->next,i++){
	os <<ItemTypes(a->type).tag <<'=';
	a->write(os);
	BAG_DEBUG(os <<" ("<<ItemTypes(a->type).typeidname <<"[" <<ItemTypes(a->type).n <<"])" <<std::endl);
	/*
	if(ItemTypes(a->type).size>0){
	  os <<'(';
	  for(i=0;i<ItemTypes(a->type).size;i++){ if(i) os <<' '; os <<a->value[i]; };
	  os <<')';
	}
	if (a->strValue.length()>0) os << " [" << a->strValue << "]";
	*/
	os <<std::endl;

      }
    }

    //! table of all item types creates so far
    static void writeExistingItemTypes(std::ostream& os){
      os <<"so far registered ItemTypes: " <<std::endl;
      for(uint i=0;i<ItemTypes.N;i++) os <<ItemTypes(i).tag <<" with type " <<ItemTypes(i).typeidname <<"[" <<ItemTypes(i).n <<"]" <<std::endl;
    }

    /*! check if there already exists a given tag with given type (typeidname & n). This 
        causes a HALT if the tag already exists but with another type */
    static int getType(const char *tag,const char *typeidname,uint n){
      uint i;
      for(i=0;i<ItemTypes.N;i++){
	if(!strcmp(ItemTypes(i).tag,tag)){
	  if(!strcmp(ItemTypes(i).typeidname,typeidname) && ItemTypes(i).n==n){
	    return i;
	  }else{
	    writeExistingItemTypes(std::cerr);
	    HALT("same tag ("<<tag<<") but different type ("<<typeidname<<"["<<n<<"]) than previously ("<<ItemTypes(i).typeidname<<"["<<ItemTypes(i).n<<"])");
	  }
	}
      }
      //that's a new type - insert it in the types table
      ItemTypes.resizeCopy(i+1);
      ItemTypes(i).tag=new char[strlen(tag)];
      strcpy(ItemTypes(i).tag,tag);
      ItemTypes(i).typeidname=typeidname;
      ItemTypes(i).n=n;
      BAG_DEBUG(std::cout <<"adding new type: " <<tag <<", " <<typeidname <<"[" <<n <<"]" <<std::endl);
      return i;
    }

    //! get the item with tag from the bag
    template<class T>
    Item* get(const char* tag,T *&x,uint n=0){
      int i;
      Item *a;
      i=getType(tag,typeid(T).name(),n);
      for(a=first;a;a=a->next) if(a->type==i) break;
      if(a) x=(T*)a->value; else x=NULL;
      return a;
    }

    //! check if the bag contains an item with a given tag and type
    bool has(const char* tag){
      Item *a;
      for(a=first;a;a=a->next) if(!strcmp(ItemTypes(a->type).tag,tag)) break;
      if(a) return true;
      return false;
    }
  };
}

//===========================================================================
//
// gnuplot calls
//

void gnuplot(const char* command);
void gnuplotEPS(const char* file,const char* command);


//===========================================================================
//
// implementations
//

#ifdef MT_IMPLEMENTATION
#  include"std.cpp"
#endif



//===========================================================================
//
// USING's
//

using std::cout;
using std::endl;
using std::ostream;
using std::istream;
using std::ofstream;
using std::ifstream;
using MT::rnd;
using MT::Color;
using MT::Parameter;
using MT::String;


#endif
