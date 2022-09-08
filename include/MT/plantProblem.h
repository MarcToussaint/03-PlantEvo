/*  Copyright (C) 2000, 2007  Marc Toussaint (mtoussai@inf.ed.ac.uk)
    under the terms of the GNU LGPL (http://www.gnu.org/copyleft/lesser.html)
    see the `std.h' file for a full copyright statement  */

#ifndef MT_plantProblem_h
#define MT_plantProblem_h

//===========================================================================
//
// the basic interface to the test function
//

/*!\brief
  the test function:
  INPUT:  string of bytes
  OUTPUT: an array of three doubles meaning (fitness, grean area, weight) as info */
void evaluatePlant(arr& info,const byteA& x);

//! the same in standard C convention (without my using my Array class)
void evaluatePlant(double info[3],byte* x,uint x_size);

/*!\brief
  display a solution:
  INPUT: string of bytes
         watch the plant interactively?
	 save a ppm image of the plant? */
void displayPlant(byteA& x,bool watch=true,char* filename=0);

//! in standard C convention (without my using my Array class)
void displayPlant(byte* x,uint x_size,bool watch=true,char* filename=0);



//===========================================================================
//
// routines for decoding a grammar
//

struct byteRule{
  byte code;
  MT::Array<byte> str;
  uint use;       //!< couter of rule applications
  uint change;    //!< a trace indicator (like eligibility traces) for mutations
  uint hierarchy; //!< the time of application during rule transformations
  uint size;
};

void readByteRule(std::istream& is, MT::Array<byteRule>& rules, arr& reals);
void writeByteRule(std::ostream& os, const MT::Array<byteRule>& rules);
void decodeRules(byteA& x,const MT::Array<byteRule>& rules,uint maxT,uint maxN);


#ifdef MT_IMPLEMENTATION
#  include"plantProblem.cpp"
#endif

#endif
