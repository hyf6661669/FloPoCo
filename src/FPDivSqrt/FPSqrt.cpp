/*
  Floating Point Square Root for FloPoCo

  Authors :
  Jeremie Detrey, Florent de Dinechin (digit-recurrence version)

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2010.
  All rights reserved.

 */



#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "FPSqrt.hpp"

using namespace std;

namespace flopoco{

#define DEBUGVHDL 0
	//#define LESS_DSPS

	FPSqrt::FPSqrt(OperatorPtr parentOp, Target* target, int wE, int wF) :
		Operator(parentOp,target), wE(wE), wF(wF) {

		ostringstream name;

		name<<"FPSqrt_"<<wE<<"_"<<wF;

		uniqueName_ = name.str();

		// -------- Parameter set up -----------------

		addFPInput ("X", wE, wF);
		addFPOutput("R", wE, wF);


		// Digit-recurrence implementation recycled from FPLibrary
		//cout << "   DDDD" <<  getTarget()->adderDelay(10) << "  " <<  getTarget()->localWireDelay() << "  " << getTarget()->lutDelay();
		vhdl << tab << declare("fracX", wF) << " <= X" << range(wF-1, 0) << "; -- fraction"  << endl;
		vhdl << tab << declare("eRn0", wE) << " <= \"0\" & X" << range(wE+wF-1, wF+1) << "; -- exponent" << endl;
		vhdl << tab << declare("xsX", 3) << " <= X"<< range(wE+wF+2, wE+wF) << "; -- exception and sign" << endl;

		vhdl << tab << declare("eRn1", wE) << " <= eRn0 + (\"00\" & " << rangeAssign(wE-3, 0, "'1'") << ") + X(" << wF << ");" << endl;

		vhdl << tab << declare(getTarget()->lutDelay(),
													 join("w",0), wF+4) << " <= \"111\" & fracX & \"0\" when X(" << wF << ") = '0' else" << endl
		     << tab << "       \"1101\" & fracX; -- pre-normalization" << endl;
		//		vhdl << tab << declare(join("d",wF+3)) << " <= '0';" << endl;
		//		vhdl << tab << declare(join("s",wF+3)) << " <= '1';" << endl;

		vhdl << tab << "-- now implementing the recurrence " << endl;
		vhdl << tab << "--  w_{i} = 2w_{i-1} -2s_{i}S_{i-1} - 2^{-i-1}s_{i}^2  for i in {1..n}" << endl;
		int maxstep=wF+2;
		for(int i=1; i<=maxstep; i++) {
		  double stageDelay= getTarget()->adderDelay(i) + 2*getTarget()->lutDelay();
			REPORT(2, "estimated delay for stage "<< i << " is " << stageDelay << "s");
			// was: int i = wF+3-step; // to have the same indices as FPLibrary
		  vhdl << tab << "-- Step " << i << endl;
		  string di = join("d", i);
		  string xi = join("x", i);
		  string wi = join("w", i);
		  string wip = join("w", i-1);
		  string si = join("s", i);
		  string sip = join("s", i-1);
		  //			string zs = join("zs", i);
		  string ds = join("ds", i);
		  string xh = join("xh", i);
		  string wh = join("wh", i);
		  vhdl << tab << declare(di) << " <= "<< wip << "("<< wF+3<<");" << endl;
		  vhdl << tab << declare(xi,wF+5) << " <= " << wip << " & \"0\";" << endl;
		  vhdl << tab << declare(ds,i+3) << " <=  \"0\" & ";
		  if (i>1)
		    vhdl 	<< sip << " & ";
		  vhdl << " (not " << di << ") & " << di << " & \"1\";" << endl;
		  vhdl << tab << declare(xh,i+3) << " <= " << xi << range(wF+4, wF+2-i) << ";" << endl;
		  vhdl << tab <<  declare(stageDelay, wh, i+3) << " <= " << xh << " - " << ds << " when " << di << "='0'" << endl
		       << tab << tab << "     else   " << xh << " + " << ds << ";" << endl;
		  vhdl << tab << declare(wi, wF+4) << " <= " << wh << range(i+1,0);
		  if(i <= wF+1)
		    vhdl << " & " << xi << range(wF+1-i, 0) << ";" << endl;
		  else
		    vhdl << ";" << endl;
		  vhdl << tab << declare(si, i) << " <= ";
		  if(i==1)
		    vhdl << "\"\" & (not " << di << ") ;"<< endl;
		  else
		    vhdl << sip /*<< range(i-1,1)*/ << " & not " << di << ";"<< endl;
		}
		string dfinal=join("d", maxstep+1);
		vhdl << tab << declare(dfinal) << " <= "<< join("w", maxstep) << of(wF+3)<<" ; -- the sign of the remainder will become the round bit" << endl;
		vhdl << tab << declare("mR", wF+3) << " <= "<< join("s", maxstep)<<" & not "<<dfinal<<"; -- result significand" << endl;

		// end of component FPSqrt_Sqrt in fplibrary
		vhdl << tab << declare(target->lutDelay(), "fR", wF+1) << " <= mR" <<range(wF, 0) << ";-- removing leading 1" << endl;
		vhdl << tab << declare("round") << " <= fR(0); -- round bit" << endl;

		vhdl << tab << declare(target->adderDelay(wF),
													 "fRn2", wF) << " <= fR" << range(wF, 1) <<" + (" << rangeAssign(wF-1, 1, "'0'") << " & round); -- rounding sqrt never changes exponents " << endl;
		vhdl << tab << declare("Rn2", wE+wF) << " <= eRn1 & fRn2;" << endl;

		vhdl << tab << "-- sign and exception processing" << endl;
		vhdl << tab <<  "with xsX select" << endl
		     << tab << tab << declare(target->lutDelay(), "xsR", 3) << " <= \"010\"  when \"010\",  -- normal case" << endl
		     << tab << tab <<  "       \"100\"  when \"100\",  -- +infty" << endl
		     << tab << tab <<  "       \"000\"  when \"000\",  -- +0" << endl
		     << tab << tab <<  "       \"001\"  when \"001\",  -- the infamous sqrt(-0)=-0" << endl
		     << tab << tab <<  "       \"110\"  when others; -- return NaN" << endl;

		vhdl << tab << "R <= xsR & Rn2; " << endl;
	}

  FPSqrt::~FPSqrt() {
  }






		void FPSqrt::emulate(TestCase * tc)
		{
			/* Get I/O values */
			mpz_class svX = tc->getInputValue("X");

			/* Compute correct value */
			FPNumber fpx(wE, wF);
			fpx = svX;
			mpfr_t x, r;
			mpfr_init2(x, 1+wF);
			mpfr_init2(r, 1+wF);
			fpx.getMPFR(x);

			if(correctRounding) {
				mpfr_sqrt(r, x, GMP_RNDN);
				FPNumber  fpr(wE, wF, r);
				/* Set outputs */
				mpz_class svr= fpr.getSignalValue();
				tc->addExpectedOutput("R", svr);
			}
			else { // faithful rounding
				mpfr_sqrt(r, x, GMP_RNDU);
				FPNumber  fpru(wE, wF, r);
				mpz_class svru = fpru.getSignalValue();
				tc->addExpectedOutput("R", svru);

				mpfr_sqrt(r, x, GMP_RNDD);
				FPNumber  fprd(wE, wF, r);
				mpz_class svrd = fprd.getSignalValue();
				/* Set outputs */
				tc->addExpectedOutput("R", svrd);
			}

			mpfr_clears(x, r, NULL);
		}




		// One test out of 4 fully random (tests NaNs etc)
		// All the remaining ones test positive numbers.
		TestCase* FPSqrt::buildRandomTestCase(int i){

			TestCase *tc;
			mpz_class a;

			tc = new TestCase(this);
			/* Fill inputs */
			if ((i & 3) == 0)
				a = getLargeRandom(wE+wF+3);
			else
				a  = getLargeRandom(wE+wF) + (mpz_class(1)<<(wE+wF+1)); // 010xxxxxx
			tc->addInput("X", a);

			/* Get correct outputs */
			emulate(tc);

			return tc;
		}

		OperatorPtr FPSqrt::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
			int wE;
			UserInterface::parseStrictlyPositiveInt(args, "wE", &wE);
			int wF;
			UserInterface::parseStrictlyPositiveInt(args, "wF", &wF);
			return new FPSqrt(parentOp, target, wE, wF);
		}


	TestList FPSqrt::unitTest(int index)
	{
		// the static list of mandatory tests
		TestList testStateList;
		vector<pair<string,string>> paramList;
		
		if(index==-1) 
		{ // The unit tests

			for(int wF=5; wF<53; wF+=1) // test various input widths
			{
					int wE = 6+(wF/10);
					while(wE>wF)
					{
						wE -= 2;
					}
					paramList.push_back(make_pair("wF",to_string(wF)));
					paramList.push_back(make_pair("wE",to_string(wE)));
					testStateList.push_back(paramList);
					paramList.clear();
			}
		}
		else     
		{
				// finite number of random test computed out of index
		}	
		return testStateList;
	}

	void FPSqrt::registerFactory(){
			UserInterface::add("FPSqrt", // name
												 "A correctly rounded floating-point square root function.",
												 "BasicFloatingPoint", // categories
												 "",
												 "wE(int): exponent size in bits; \
wF(int): mantissa size in bits",
												 "",
												 FPSqrt::parseArguments,
												 FPSqrt::unitTest
												 ) ;

		}
	}
