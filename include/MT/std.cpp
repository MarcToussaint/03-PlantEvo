/*  Copyright (C) 2000, 2006  Marc Toussaint (mtoussai@inf.ed.ac.uk)
    under the terms of the GNU LGPL (http://www.gnu.org/copyleft/lesser.html)
    see the `std.h' file for a full copyright statement  */

#include "std.h"

#include <unistd.h>

#ifdef MT_QT
#  include<qapplication.h>
#endif

#ifndef MT_ConfigFileName
#  define MT_ConfigFileName "MT.cfg"
#  define MT_LogFileName "MT.log"
#endif


//===========================================================================
//
// Bag container
//

MT::Mem<MT::Bag::ItemType> MT::Bag::ItemTypes;

char* MT::String::readSkipSymbols=" \t";
char* MT::String::readStopSymbols="\n\r";
bool  MT::String::eatStopSymbol  = false;


//===========================================================================
//
// other implementations
//




//===========================================================================
//
// utilities in MT namespace
//

namespace MT{
  int argc;
  char** argv;
  std::ifstream cfgFile;
  bool cfgOpenFlag=false;
  bool cfgLock=false;
  Mem<void*> parameters;
  bool IOraw(false);
  bool noLog(true);

#ifndef MT_TIMEB
  timeval startTime;
#else
  _timeb startTime;
#endif
  double timerTime;

#ifdef MT_QT
  bool qAppCreatedFlag=false;
#endif

  class Demon{
  public:
    std::ofstream logFile;
    int logstat;

    Demon(){
      logstat=0;
      timerTime=MT::realTime();
#ifndef MT_TIMEB
      gettimeofday(&startTime,0);
#else
      _ftime(&startTime);
#endif
    }
    ~Demon(){
      if(logstat){ //don't open a log file anymore in the destructor
	MT::log() <<"Execution stop:      " <<date()
	  <<"Ellapsed real time:  " <<realTime() <<"sec\n"
	  <<"Process  user time:   "<<cpuTime() <<"sec" <<std::endl;
#ifndef MT_TIMEB
	MT::log() <<"Process system time: " <<sysTime() <<"sec" <<std::endl;
#endif
	MT::log().close();
      }
#ifdef MT_QT
      if(qAppCreatedFlag){
	qApp->processEvents();
	qApp->quit();
      }
#endif
    }

    std::ofstream& log(const char *name){
      if(!logstat && noLog) logstat=1;
      if(!logstat){
	logFile.open(name);
	if(!logFile.good()) MT_MSG("could not open log-file `"<<name<<"' for output");
	logstat=1;
      }
      return logFile;
    }
  }demon;

  std::ofstream& log(const char *name){ return demon.log(name); }

  bool contains(const char* s,char c){
    for(uint i=0;s[i];i++) if(s[i]==c) return true;
    return false;
  }

  char skip(std::istream& is,char* skipchars){
    char c;
    do{ c=is.get(); } while(contains(skipchars,c));
    is.putback(c);
    return c;
  }

  void skipLine(std::istream& is){
    char c;
    do{ c=is.get(); } while(c!='\n');
  }

  byte bit(byte *str,uint i){ return (str[i>>3] >> (7-(i&7))) & 1; }
  //byte bit(byte b,uint i){ return (b >> (7-(i&7))) & 1; }
  //void set(byte *state,uint i){ state[i>>3] |= 1 << (7-(i&7)); }
  //void del(byte *state,uint i){ state[i>>3] &= (byte)~(1 << (7-(i&7))); }
  //void flip(byte *str,uint i){ str[i>>3] ^= 1 << (7-(i&7)); }

  void flip(byte& b,uint i){ b ^= 1 << (7-(i&7)); }

  void flip(int& b,uint i){ b ^= 1 << (7-(i&7)); }

  double modMetric(double x,double y,double mod){
    double d=fabs(x-y);
    d=fmod(d,mod);
    if(d>mod/2.) d=mod-d;
    return d;
  }

  double sign(double x){ if(x<0.) return -1.; return 1.; }

  double linsig(double x){ if(x<0.) return 0.; if(x>1.) return 1.; return x; }

  void constrain(double& x,double a,double b){ if(x<a) x=a; if(x>b) x=b; }

  double phi(double x,double y){
    if(x==0. || ::fabs(x)<1e-10) if(y>0.) return MT_PI/2.; else return -MT_PI/2.;
    double p=::atan(y/x);
    if(x<0.) if(y<0.) p-=MT_PI; else p+=MT_PI;
    if(p>MT_PI)  p-=2.*MT_PI;
    if(p<-MT_PI) p+=2.*MT_PI;
    return p;
  }
  double dphi(double x,double y,double dx,double dy){
    //return (dy*x - dx*y)/sqrt(x*x+y*y);
    if(x==0. || ::fabs(x)<1e-10) if(y>0.) return -dx/y; else return dx/y;
    double f=y/x;
    return 1./(1.+f*f)*(dy/x - f/x*dx);
  }

  double DIV(double x,double y,bool force){
    if(x==0.) return 0.;
    if(force){ if(y==0.) return 0.; }else CHECK(y!=0,"Division by Zero!");
    return x/y;
  }

  double Log(double x){
    if(x<.001) x=.001; return ::log(x);
  }

  uint Log2(uint n){
    uint l=0;
    n=n>>1;
    while(n){ l++; n=n>>1; }
    return l;
  }

  double sqr(double x){ return x*x; }

  double sinc(double x){
    if(fabs(x)<1e-10) return 1.-.167*x*x;
    return ::sin(x)/x;
  }

  double cosc(double x){
    if(fabs(x)<1e-10) return 1.-.167*x*x;
    return ::cos(x)/x;
  }

#define AXETS 1280
#define AXETR 10.
  double approxExp(double x){
    static bool initialized=false;
    static double ExpTable [AXETS]; //table ranges from x=-10 to x=10
    int i;
    if(!initialized){
      for(i=0;i<AXETS;i++) ExpTable[i]=::exp(AXETR*(2*i-AXETS)/AXETS);
      initialized=true;
    }
    x*=.5*AXETS/AXETR;
    i=(int)x;
    x-=(double)i; //x = residual
    i+=AXETS/2; //zero offset
    if(i>=AXETS-1) return ExpTable[AXETS-1];
    if(i<=0) return 0.;//ExpTable[0];
    return (1.-x)*ExpTable[i] + x*ExpTable[i+1];
  }

  void open(std::ofstream& fs,const char *name,char* errmsg){
    fs.clear();
    fs.open(name);
    log() <<"opening output file `" <<name <<"'" <<std::endl;
    if(!fs.good()) MT_MSG("could not open file `" <<name <<"' for output"<<errmsg);
  }
  void open(std::ifstream& fs,const char *name,char* errmsg){ 
    fs.clear();
    fs.open(name);
    log() <<"opening input file `" <<name <<"'" <<std::endl;
    if(!fs.good()) HALT("could not open file `" <<name <<"' for input"<<errmsg);
  }
  double realTime(){
#ifndef MT_TIMEB
    static timeval t; gettimeofday(&t,0);
    return ((double)(t.tv_sec-startTime.tv_sec-1) +
	    (double)((long)1000000+t.tv_usec-startTime.tv_usec)/1000000.);
#else
    static _timeb t; _ftime(&t);
    return ((double)(t.time-startTime.time-1) +
	    (double)((unsigned short)1000+t.millitm-startTime.millitm)/1000.);
#endif
  }

  double cpuTime(){
#ifndef MT_TIMEB
    static tms t; times(&t);
    return ((double)t.tms_utime)/sysconf(_SC_CLK_TCK);
#else
    return ((double)clock())/CLOCKS_PER_SEC; //MSVC: CLOCKS_PER_SEC=1000
#endif
  }

  double sysTime(){
#ifndef MT_TIMEB
    static tms t; times(&t);
    return ((double)(t.tms_stime))/sysconf(_SC_CLK_TCK);
#else
    HALT("sysTime() is not implemented for Windows!");
    return 0.;
#endif
  }

  double totalTime(){
#ifndef MT_TIMEB
    static tms t; times(&t);
    return ((double)(t.tms_utime+t.tms_stime))/sysconf(_SC_CLK_TCK);
#else
    HALT("totalTime() is not implemented for Windows!");
    return 0.;
#endif
  }

  char *date(){ static time_t t; time(&t); return ctime(&t); }

  void resetTimer(){ timerTime=realTime(); }

  double getTimer(bool reset,double relTime){
    if(relTime==-1.) relTime=timerTime;
    double c=realTime()-relTime;
    log() <<"Reading timer after user REAL time: " <<c <<"sec\n";
    if(reset) timerTime=realTime(); 
    return c;
  }

  //void wait(double sec){ double t=realTime(); while(realTime()-t<sec); }
  void wait(double sec){
#ifndef MT_TIMEB
    timeval tv;
    tv.tv_sec = (int)(floor(sec));
    sec -= (double)tv.tv_sec;
    tv.tv_usec = (int)(floor(1000000. * sec));
    int r = select(1, NULL, NULL, NULL, &tv);
    if(r==-1) HALT("select() failed");
    /* r=0 time is up
       r!=0 data in NULL stream available (nonsense)
       */
#else
#  ifdef _WINDOWS_
    MsgWaitForMultipleObjects( 0, NULL, FALSE, (int)(1000.*sec), QS_ALLEVENTS);
#  else
    double t=realTime(); while(realTime()-t<sec);
#  endif
#endif
  }

  void wait(){
    char c;
    std::cout <<"<hit key>" <<std::endl;
    std::cin.getline(&c,1);
  }

  long mem(){
#ifndef MT_TIMEB
    static rusage r; getrusage(RUSAGE_SELF,&r);
    return r.ru_idrss;
#else
    HALT("MT::mem() is not implemented for Windows!");
    return 0;
#endif
  }

  void init(int _argc, char *_argv[]){
    argc=_argc; argv=_argv;
    time_t t; time(&t);
    char *name;
    if(checkOption("nolog")) noLog=true; else noLog=false;
    if(!getOption("log",name)) name=MT_LogFileName;
    log(name);
    if(!log().good()) MT_MSG(" -- use `-nolog' or `-log' option to specify the log file");
    log() <<"Compiled at:     " <<__DATE__ <<" "<<__TIME__ <<"\n";
    log() <<"Execution start: " <<ctime(&t);
    log() <<"Program call:    '"; for(int i=0;i<argc;i++) log() <<argv[i] <<" ";
    log() <<"\b'" <<std::endl;
  }

  bool checkOption(const char *tag){
    for(int n=1;n<argc;n++) if(argv[n][0]=='-' && !strcmp(tag,argv[n]+1)){
      return true;
    }
    return false;
  }

  bool getOption(const char *tag,char* &option){
    for(int n=1;n+1<argc;n++) if(argv[n][0]=='-' && !strcmp(tag,argv[n]+1)){
      option=argv[n+1];
      return true;
    }
    return false;
  }

  void openConfigFile(char* name){
    log() <<"opening config file ";
    if(!name) if(!getOption("cfg",name)) name=MT_ConfigFileName;
    if(cfgOpenFlag){ cfgFile.close(); log() <<"(old config file closed) "; }
    log() <<"'" <<name <<"'";
    cfgFile.clear();
    cfgFile.open(name);
    cfgOpenFlag=true;
    if(!cfgFile.good()){
      //MT_MSG("couldn't open config file "<<name);
      log() <<" - failed";
    }
    log() <<std::endl;
  }

#ifdef MT_QT
  InitQt::InitQt(){
    if(!qAppCreatedFlag){
      int ARGN=2;
      char *ARGV[2]={"x","-style=windows"};
      new QApplication(ARGN,ARGV);
      qAppCreatedFlag=true;
    }
  }
#endif
}


//===========================================================================
//
// random number generator
//

namespace MT{ Rnd rnd; }

unsigned long MT::Rnd::seed(unsigned long n){
  unsigned long s,c;
  if(n>12345){ s=n; c=n%113; }
  else{ s=12345; c=n; }
  while(c--) s*=65539;
  s=s>>1;
  seed250(s);
  ready=true;
  return n;
}

unsigned long MT::Rnd::seed(){
  Parameter<uint> _seed("seed",0);
  return seed(_seed);
}

unsigned long MT::Rnd::clockSeed(){
  unsigned long s;
#ifndef MT_TIMEB
  timeval t; gettimeofday(&t,0); s=1000000L*t.tv_sec+t.tv_usec;
#else
  _timeb t; _ftime(&t); s=1000L*t.time+t.millitm;
#endif
  log() <<"random clock seed: " <<s <<std::endl;
  return seed(s);
}

double MT::Rnd::gauss(){
  double w,v,rsq,fac;
  do{
    v   = 2 * uni() - 1;
    w   = 2 * uni() - 1;
    rsq = v*v + w*w;
  } while( rsq >= 1 || rsq == 0 );
  fac  = ::sqrt( -2 * ::log( rsq ) / rsq );
  return v*fac;
}

uint MT::Rnd::poisson(double mean){
  if(mean>100){
    int i=(int)::floor(mean+::sqrt(mean)*gauss()+.5);
    return (i>0)?(uint)i:0;
  }
  uint count = 0;
  double bound, product;
  if(mean>=0){
    bound=::exp(-mean);
    product=uni();
    while(product>=bound){
      count++;
      product*=uni( );
    }
  }
  return count;
}

void  MT::Rnd::seed250(long seed){
  int      i; 
  long     j, k;
  
  if(seed<=0) seed=1;
    
  for(i=0;i<250;++i){           // Schleife ueber Zufallsfeld
    k = seed / 127773;          // Modulozufallszahlengenerator
    seed = 16807 * ( seed - k*127773 ) - 2836 * k;
    if(seed<0) seed += 0x7FFFFFFF;
    rfield[i] = seed;
  }

  // Masken ueberlagern
  k = 0x7FFFFFFF;       
  j = 0x40000000;
  for(i=1;i<250;i+=8) rfield[i] = (rfield[i] & k) | j;
    
  // rpoint initialisieren
  rpoint = 249;            
    
  // Anfangszahlen verwerfen
  for(i=0;i<4711;++i) rnd250();
}    

void gnuplot(const char* command){
  std::ofstream os("z.plt");
#ifndef MT_MSVC
  os
    <<"set terminal x11\n" //persist noraise\n"
    <<"set data style lines\n"
    <<command <<std::endl
    <<"pause mouse" <<std::endl;
  os.close();
  if(system("gnuplot z.plt"))
    MT_MSG("attempt to call gnuplot on via `gnuplot z.plt' failed");
#else
  os
    <<"set data style lines\n"
    <<command <<std::endl
    <<"pause mouse" <<std::endl;
  os.close();
  if(system("wgnuplot z.plt"))
    MT_MSG("attempt to call gnuplot on windows via `wgnuplot z.plt' failed (probably not installed)");
#endif
}

void gnuplotEPS(const char* filename,const char* command){
  std::ofstream os("z.plt");
  os
    <<"set terminal postscript 24 dl 5; set key samplen 2\n"
    <<"set output '" <<filename <<"'\n"
    <<"set data style lines\n"
    <<command <<std::endl;
  os.close();
#ifndef MT_MSVC
  if(system("gnuplot z.plt"))
    MT_MSG("attempt to call gnuplot on via `gnuplot z.plt' failed - given command:\n"<<command);
#else
  if(system("wgnuplot z.plt"))
    MT_MSG("attempt to call gnuplot on windows via `wgnuplot z.plt' failed (probably not installed) - given command:\n"<<command);
#endif
}

//===========================================================================
//
// Cumulative probability for the Standard Normal Distribution 
//

namespace MT{
  double erf(double x){
    double t, z, retval;
    z = fabs( x );
    t = 1.0 / ( 1.0 + 0.5 * z );
    retval = t * exp( -z * z - 1.26551223 + t *
		    ( 1.00002368 + t *
		    ( 0.37409196 + t *
		    ( 0.09678418 + t *
		    ( -0.18628806 + t *
		    ( 0.27886807 + t *
		    ( -1.13520398 + t *
		    ( 1.48851587 + t *
		    ( -0.82215223 + t *
		    0.1708727 ) ) ) ) ) ) ) ) );
    if( x < 0.0 ) return retval - 1.0;
    return 1.0 - retval;
  }
  
  //! the integral of N(0,1) from -infty to x
  double gaussInt(double x){
    return .5*(1.+erf(x/MT_SQRT2));
  }

  //! expectation \f$\int_x^\infty {\cal N}(x) x dx\f$ when integrated from -infty to x
  double gaussIntExpectation(double x){
    double norm=gaussInt(x) / (::sqrt(MT_2PI));
    return - norm*MT::approxExp(-.5*x*x);
  }
}



