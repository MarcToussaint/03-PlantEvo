/*  Copyright (C) 2000, 2006  Marc Toussaint (mtoussai@inf.ed.ac.uk)
    under the terms of the GNU LGPL (http://www.gnu.org/copyleft/lesser.html)
    see the `std.h' file for a full copyright statement  */

#include "list.h"

void random(doubleAL& L,uint N,uint d,double lo,double hi){
  uint i;
  L.clear();
  doubleA a;
  for(i=0;i<N;i++){
    a.resize(d);
    rndUni(a,lo,hi);
    L.append(a);
  }
  L.weightEqual();
}

void random(intAL& L,uint N,uint d,int lo,int hi){
  uint i;
  L.clear();
  intA a;
  for(i=0;i<N;i++){
    a.resize(d);
    rndInt(a,lo,hi);
    L.append(a);
  }
  L.weightEqual();
}
