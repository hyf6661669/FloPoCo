/*
  Floating Point Normal CDF for FloPoCo
 
  Authors : 
  David Thomas

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  Imperial
  2008-2013.
  All rights reserved.

 */
 


#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <gmp.h>
#include <mpfr.h>

#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "IntAdder.hpp"
#include "IntMultiplier.hpp"
#include "IntSquarer.hpp"
#include "FPSquarer.hpp"
#include "FPExp.hpp"
#include "FPMultiplier.hpp"
#include "FP2Fix.hpp"
#include "FP2FixV2.hpp"
#include "Fix2FP.hpp"
#include "FixFunctions/FunctionEvaluator.hpp"

#include "FPNormalCDF.hpp"

#include "UtilSollya.hh"

using namespace std;

namespace flopoco{

	void FPNormalCDF::NCD(mpfr_t r, mpfr_t x, mpfr_rnd_t mode)
	{
		mpfr_t w;
		mpfr_init2(w, getToolPrecision());
		
		mpfr_set_ui(w, 2u, MPFR_RNDN);		// w=2
		mpfr_sqrt(w,w,MPFR_RNDN);								// w=sqrt(2)
		mpfr_div(w, x, w, MPFR_RNDN);							// w=x/sqrt(2)
		mpfr_neg(w, w, MPFR_RNDN);							// w=-x/sqrt(2)
		mpfr_erfc(w, w, MPFR_RNDN);							// w=erfc(-x/sqrt(2))
		mpfr_div_ui(w, w,2u, MPFR_RNDN);					// w=rfc(-x/sqrt(2))/2
		
		mpfr_set(r, w, mode);
		
		mpfr_clear(w);
	}
	
	void FPNormalCDF::NCDInv(mpfr_t r, mpfr_t p)
	{
		mpfr_t lowerX, upperX, lowerP, upperP, midX, midP;
		mpfr_inits2(getToolPrecision(), lowerX, upperX, lowerP, upperP, midX, midP, NULL);
		
		int tries=0;
		mpfr_set_si(lowerX, -1, MPFR_RNDN);
		do{
			mpfr_mul_si(lowerX, lowerX, 2, MPFR_RNDN);
			NCD(lowerP, lowerX);
			//REPORT(FULL, "NCDInv : lowerX="<<mpfr_get_d(lowerX,MPFR_RNDN)<<" -> lowerP = "<<mpfr_get_d(lowerP,MPFR_RNDN));
			if(++tries>20){
				exit(1);
			}
		}while(mpfr_cmp(p,lowerP)<=0);
		
		mpfr_set_si(upperX, +1, MPFR_RNDN);
		do{
			mpfr_mul_si(upperX, upperX, 2, MPFR_RNDN);
			NCD(upperP, upperX);
			//REPORT(FULL, "NCDInv : upperX="<<mpfr_get_d(upperX,MPFR_RNDN)<<" -> upperP = "<<mpfr_get_d(upperP,MPFR_RNDN));
		}while(mpfr_cmp(p,upperP)>0);
		
		do{
			assert(mpfr_cmp(lowerX,upperX)<0);
			assert(mpfr_cmp(lowerP,p)<0);
			assert(mpfr_cmp(p,upperP)<=0);
			
			mpfr_add(midX, lowerX, upperX, MPFR_RNDU);
			mpfr_div_si(midX, midX, 2, MPFR_RNDN);
			
			//REPORT(FULL, "  lowerX="<<mpfr_get_d(lowerX, MPFR_RNDN)<<", upperX="<<mpfr_get_d(upperX, MPFR_RNDN));
			//REPORT(FULL, "  lowerP="<<mpfr_get_d(lowerP, MPFR_RNDN)<<", upperP="<<mpfr_get_d(upperP, MPFR_RNDN));
			
			if(mpfr_cmp(midX,upperX)==0){
				break;
			}
			
			NCD(midP, midX);
			
			if(mpfr_cmp(midP,p)<0){
				mpfr_set(lowerP, midP,MPFR_RNDN);
				mpfr_set(lowerX, midX, MPFR_RNDN);
			}else{
				mpfr_set(upperP, midP, MPFR_RNDN);
				mpfr_set(upperX, midX, MPFR_RNDN);
			}
		}while(1);
		
		REPORT(FULL, "NCDInv("<<mpfr_get_d(p,MPFR_RNDN)<<") = "<<mpfr_get_d(upperX,MPFR_RNDN));
		
		mpfr_set(r, upperX, MPFR_RNDN);
		mpfr_clears(lowerX, midX, upperX, lowerP, midP, upperP, NULL);
	}
	
	void FPNormalCDF::B(mpfr_t r, mpfr_t x, mpfr_rnd_t mode)
	{
		mpfr_t w;
		
		mpfr_init2(w, getToolPrecision());
		mpfr_sqr(w, x, MPFR_RNDN);	// w=x^2
		mpfr_div_si(w, w, -2, MPFR_RNDN);	// w=-x^2/2
		mpfr_exp(r, w, mode);	// r=exp(-x^2/2)
		
		mpfr_clear(w);
	}
	
	void FPNormalCDF::F(mpfr_t r, mpfr_t x, mpfr_rnd_t mode)
	{
		mpfr_t NCDx, Bx;
		
		mpfr_init2(NCDx, getToolPrecision());
		NCD(NCDx, x);
		
		mpfr_init2(Bx, getToolPrecision());
		B(Bx, x);
		
		mpfr_div(r, NCDx, Bx, mode);
		
		mpfr_clear(NCDx);
		mpfr_clear(Bx);		
	}
	
	// Calculate the minimum input value we have to worry about. Below this point
	// it will underflow to zero
	void FPNormalCDF::CalcMinX()
	{
		// First look at the minimum by input range.
		FPNumber minInput(wInE, wInF, FPNumber::largestNegative);
		
		// Then see if the users range will help
		mpfr_t minInputX;
		mpfr_init2(minInputX, wInF);
		minInput.getMPFR(minInputX);
		REPORT(DETAILED, "minInput from input range="<<mpfr_get_d(minInputX, MPFR_RNDN));
		
		if(mpfr_cmp_d(minInputX, minInputValue)<0){
			mpfr_set_d(minInputX, minInputValue, MPFR_RNDN);
			REPORT(DETAILED, "minInput updated from user="<<mpfr_get_d(minInputX, MPFR_RNDN));	
		}
		
		// Now work backwards from output representation
		FPNumber minProb(wOutE, wOutF, FPNumber::smallestPositive);
		mpfr_t minProbP, minProbX;
		mpfr_init2(minProbP, getToolPrecision());
		mpfr_init2(minProbX, getToolPrecision());
		minProb.getMPFR(minProbP);
		REPORT(DEBUG, "Working out input range from min probability of "<<mpfr_get_d(minProbP, MPFR_RNDN));
		NCDInv(minProbX, minProbP);
		mpfr_prec_round(minProbX, wInF, MPFR_RNDN);
		if(mpfr_cmp(minInputX, minProbX)<0){
			mpfr_set(minInputX, minProbX, MPFR_RNDN);
			REPORT(DETAILED, "minInput updated from output probability="<<mpfr_get_d(minInputX, MPFR_RNDN));
		}

		mpfr_init2(m_minInputX, wInF);
		mpfr_set(m_minInputX, minInputX, MPFR_RNDN);
		
		mpfr_clears(minInputX, minProbX, minProbP, NULL);
	}
	
	int FPNormalCDF::CalcFxOutputLSB(int relativeBitsNeeded)
	{
		mpfr_t tmp;
		mpfr_init2(tmp, getToolPrecision());
		
		mpfr_set(tmp, m_minInputX, MPFR_RNDN);
		F(tmp, tmp);
		
		REPORT(DEBUG, "  F("<<-mpfr_get_d(m_minInputX,MPFR_RNDN)<<") = "<<mpfr_get_d(tmp, MPFR_RNDN));
		
		mpfr_log2(tmp, tmp, MPFR_RNDN);
		mpfr_floor(tmp, tmp);
		
		REPORT(DEBUG, "  floor(log2(ans))="<<mpfr_get_d(tmp, MPFR_RNDN));
		
		int res=mpfr_get_si(tmp, MPFR_RNDN);
		
		mpfr_clear(tmp);
		
		return res-relativeBitsNeeded;
	}
	
	int FPNormalCDF::CalcFxInputMSB()
	{
		mpfr_t tmp;
		mpfr_init2(tmp, getToolPrecision());
		
		mpfr_si_sub(tmp, 0, m_minInputX, MPFR_RNDN);
		mpfr_log2(tmp, tmp, MPFR_RNDN);
		mpfr_ceil(tmp, tmp);
		
		int res=mpfr_get_si(tmp, MPFR_RNDN);
		mpfr_clear(tmp);
		return res;		
	}

	FPNormalCDF::FPNormalCDF(Target* target, int wInE, int wInF, int wOutE, int wOutF, double minInputValue) :
		Operator(target),
		wInE(wInE), wInF(wInF),
		wOutE(wOutE), wOutF(wOutF),
		minInputValue(minInputValue),
		funcFx(NULL)
	{
		m_debugOutputs=true;

		ostringstream name;

		name<<"FPNormalCDF_"<<wInE<<"_"<<wInF<<"_"<<wOutE<<"_"<<wOutF<<"_n"<<int(-minInputValue);
		setName(name.str()); 
		setCopyrightString("David B. Thomas (2013)");
		
		addFPInput ("X", wInE, wInF);
		addFPOutput("R", wOutE, wOutF);
		
		REPORT(INFO, "Creating FPNormalCDF with input=("<<wInE<<","<<wInF<<"), output=("<<wOutE<<","<<wOutF<<"), minInputValue="<<minInputValue<<"\n");
		
		CalcMinX();
		
		/* We're going to use:
			NCD(x) = normcdf(x)
			B(x) = exp(-x^2/2)
			F(x) = NCD(x) / B(X)
		So:
			NCDa(x) = F(x) * B(x)
		*/		
		
		
		
		////////////////////////////////////////////////////////////////
		// Part 1 : Calculate Bx=B(x)
		REPORT(DETAILED, "Building up B(x) part.");
		
		// Do sqr(X)=X^2
		wSqrXF=wOutF+8;
		int wSqrX=2+1+wInE+wSqrXF;
		FPSquarer *opSquarer=new FPSquarer(target, wInE, wInF, wSqrXF);
		oplist.push_back(opSquarer);
		inPortMap(opSquarer, "X", "X");
		outPortMap(opSquarer, "R", "sqrX");
		vhdl<<tab<<instance(opSquarer, "squarer")<<"\n";
		syncCycleFromSignal("sqrX");
		
		// negHalfSqrX = -x^2/2
		vhdl<<tab<<declare("sqrX_flags", 2)<<" <= sqrX"<<range(wSqrX-1,wSqrX-2)<<";\n";
		vhdl<<tab<<declare("sqrX_expnt", wInE)<<" <= sqrX"<<range(wInE+wSqrXF-1,wSqrXF)<<";\n";
		vhdl<<tab<<declare("sqrX_frac", wSqrXF)<<" <= sqrX"<<range(wSqrXF-1,0)<<";\n";
		// TODO : Handle under-flow of numbers
		vhdl<<tab<<declare("negHalfSqrX",wSqrX)<<" <= sqrX when not sqrX_flags=\"01\" else sqrX_flags & '1' & (sqrX_expnt-1) & sqrX_frac;\n";
		syncCycleFromSignal("negHalfSqrX");
		
		// Bx=exp(-abs(x)^2/2)
		wBxF=wSqrXF;
		wBxE=wInE;
		int wB=2+1+wBxE+wBxF;
		FPExp *opExp=new FPExp(target, wInE, wSqrXF, 0, 0);
		oplist.push_back(opExp);
		inPortMap(opExp, "X", "negHalfSqrX");
		outPortMap(opExp, "R", "Bx");
		vhdl<<tab<<instance(opExp, "exp")<<"\n";
		syncCycleFromSignal("Bx");
		
		//////////////////////////////////////////////////////
		// Part 2 : Calculate Fx=F(x)
		REPORT(DETAILED, "Building up F(x) part.");
		
		setCycleFromSignal("X");
		
		/* First produce absX=abs(X), and reflect_at_end=X<0 */
		vhdl<<tab<<declare("reflect_at_end")<<" <= X("<<wInE+wInF<<");\n";
		vhdl<<tab<<declare("absX",2+1+wInE+wInF)<<" <= X"<<range(2+1+wInE+wInF-1,1+wInE+wInF)<<"&'0'&X"<<range(wInE+wInF-1,0)<<";\n";
		syncCycleFromSignal("absX");
		
		// Convert it to fixed-point.
		fixXLSB=-(wInF+8);
		fixXMSB=CalcFxInputMSB()+1;
		REPORT(INFO, "fixX: MSB=2^"<<fixXMSB<<", fixXLSB=2^"<<fixXLSB);
		// TODO : For some reason FP2Fix wants us to have more fractional bits in the input than there are in
		// the output. I have no idea why this should be the case...
		// TODO : Actually, it needs one _more_ bit as well.
		// TODO: and another one! But I think that is my fault
		int extraBits=0; //(-fixXLSB)-wInF+1+1;
		REPORT(DETAILED, "wInF="<<wInF<<", -fixXLSB="<<-fixXLSB);
		if(extraBits<=0){
			vhdl<<tab<<declare("absX_Hack", 3+wInE+wInF)<<" <= absX;\n";
		}else{
			vhdl<<tab<<declare("absX_Hack", 3+wInE+wInF+extraBits)<<" <= absX & "<<zg(extraBits)<<";\n";
			REPORT(DETAILED, "  adding "<<extraBits<<" to keep FP2Fix happy");
		}
		syncCycleFromSignal("absX_Hack");
		
		// TODO: There is something very fishy about this whole conversion, I think I am wasting LSBs here
		FP2FixV2 *opFP2Fix=new FP2FixV2(target, /*LSB0*/ fixXLSB-1, /*MSBO*/ fixXMSB-1, /*Signed*/ false, /*wER*/ wInE, wInF+extraBits, /*trunc_p*/ false);
		oplist.push_back(opFP2Fix);
		inPortMap(opFP2Fix, "I", "absX_Hack");
		outPortMap(opFP2Fix, "O", "fixX");
		vhdl<<tab<<instance(opFP2Fix, "fp2fix")<<"\n";
		syncCycleFromSignal("fixX");
		
		// Do the actual evaluation of F(x)
		FxLSB=CalcFxOutputLSB(wOutF+2);
		FxMSB=-1;	// The maximum is always 0.5
		int degree=4;	// TODO
		std::string funcStr;
		{
			std::stringstream xStr;
			xStr<<"(x*2^"<<fixXMSB<<")";
			std::string x=xStr.str();
			std::stringstream fStr;
			fStr<<" erfc("<<x<<"/sqrt(2))/2/exp(-"<<x<<"^2/2),0,1,1";
			funcStr=fStr.str();
			
			funcFx=new Function(funcStr);
			
			/*for(double i=0;i<1.0;i+=0.125){
				REPORT(DEBUG, " f("<<i<<") -> "<<funcFx->eval(i));
			}*/
		}
		REPORT(INFO, "funcStr="<<funcStr);
		REPORT(INFO, "  Out, lsb=2^"<<FxLSB);
		opFx=new FunctionEvaluator(target, funcStr.c_str(), fixXMSB-fixXLSB+1, -FxLSB, degree);
		oplist.push_back(opFx);
		inPortMap(opFx, "X", "fixX");
		outPortMap(opFx, "R", "fixFx");
		vhdl<<tab<<instance(opFx, "FxEval")<<"\n";
		syncCycleFromSignal("fixFx");
		
		REPORT(INFO, "Finished building function evaluator.");
		
		// Convert it back to floating-point
		wFxF=wOutF+2;
		wFxE=wOutE;
		// The function evaluator may have produced more bits than we want
		vhdl<<tab<<declare("fixFx_trunc", FxMSB-FxLSB+1)<<" <= fixFx"<<range(FxMSB-FxLSB,0)<<";\n";
		
		Fix2FP *opFix2FP=new Fix2FP(target, FxLSB, FxMSB+1, /*Signed*/ false, wFxE, wFxF);
		
		oplist.push_back(opFix2FP);
		inPortMap(opFix2FP, "I", "fixFx_trunc");
		outPortMap(opFix2FP, "O", "Fx");
		vhdl<<tab<<instance(opFix2FP, "fix2fp")<<"\n";
		
		///////////////////////////////////////////////////
		// Part 3: Join them back up and handle any special cases
		REPORT(DETAILED, "Joining two sides.");
		
		syncCycleFromSignal("Fx");
		syncCycleFromSignal("Bx");
		nextCycle();
		
		vhdl<<tab<<declare("Fx_in",3+wFxE+wFxF)<<" <= Fx;\n";
		vhdl<<tab<<declare("Bx_in",3+wBxE+wBxF)<<" <= Bx;\n";
		
		FPMultiplier *opMult=new FPMultiplier(target, wBxE, wBxF, wFxE, wFxF, wOutE, wOutF);
		
		oplist.push_back(opMult);
		inPortMap(opMult, "X", "Bx_in");
		inPortMap(opMult, "Y", "Fx_in");
		outPortMap(opMult, "R", "result");
		vhdl<<tab<<instance(opMult, "finalMult")<<"\n";
		
		
		syncCycleFromSignal("result");
		
		if(m_debugOutputs){
			addFPOutput("oDebug_sqrX", wInE, wSqrXF);
			vhdl<<tab<<"oDebug_sqrX <= sqrX;\n";
			//addOutput("oDebug_sqrX_frac", wSqrXF);
			//vhdl<<tab<<"oDebug_sqrX_frac <= sqrX"<<range(wSqrXF-1,0)<<";\n";
			//addOutput("oDebug_sqrX_expnt", wInE);
			//vhdl<<tab<<"oDebug_sqrX_expnt <= sqrX"<<range(wInE+wSqrXF-1,wSqrXF)<<";\n";
			
			addFPOutput("oDebug_negHalfSqrX", wInE, wSqrXF);
			vhdl<<tab<<"oDebug_negHalfSqrX <= negHalfSqrX;\n";
			//addOutput("oDebug_negHalfSqrX_frac", wSqrXF);
			//vhdl<<tab<<"oDebug_negHalfSqrX_frac <= negHalfSqrX"<<range(wSqrXF-1,0)<<";\n";
			//addOutput("oDebug_negHalfSqrX_expnt", wInE);
			//vhdl<<tab<<"oDebug_negHalfSqrX_expnt <= negHalfSqrX"<<range(wInE+wSqrXF-1,wSqrXF)<<";\n";
			
			addFPOutput("oDebug_Bx", wBxE, wBxF);
			vhdl<<tab<<"oDebug_Bx <= Bx_in;\n";
			//addOutput("oDebug_Bx_expnt", wBxE);
			//vhdl<<tab<<"oDebug_Bx_expnt <= Bx"<<range(wBxE+wBxF-1,wBxF)<<";\n";
			//addOutput("oDebug_Bx_frac", wBxF);
			//vhdl<<tab<<"oDebug_Bx_frac <= Bx"<<range(wBxF-1,0)<<";\n";
			
			addOutput("oDebug_fixX", fixXMSB-fixXLSB+1);
			vhdl<<tab<<"oDebug_fixX <= fixX;\n";
			
			addOutput("oDebug_fixFx_trunc", FxMSB-FxLSB+1);
			vhdl<<tab<<"oDebug_fixFx_trunc <= fixFx_trunc;\n";
			
			addFPOutput("oDebug_Fx", wFxE,wFxF);
			vhdl<<tab<<"oDebug_Fx <= Fx_in;\n";
		}
		
		vhdl<<tab<<"R <= result;\n";
		
		REPORT(DETAILED, "Done.");
	}
  
  FPNormalCDF::~FPNormalCDF() {
	  mpfr_clear(m_minInputX);
	  if(funcFx)
		delete funcFx;
  }
  
	void FPNormalCDF::buildStandardTestCases(TestCaseList* tcl)
	{
		{
			mpfr_t a,b;
			mpfr_init2(a, 1+wOutF);
			mpfr_init2(b, 1+wOutF);
			
			gmp_randstate_t rng;
			gmp_randinit_default(rng);
			
			for(unsigned i=0;i<1000;i++){
				mpfr_urandomb(a, rng);
				FPNumber fpr(wOutE, wOutF, a);
				fpr.getMPFR(b);
				if(mpfr_cmp(a,b)){
					mpfr_printf("Had %bR -> %bR\n", a, b);
					exit(1);
				}
			}
			
			gmp_randclear(rng);
		}
		
		mpfr_t x;
		
		mpfr_init2(x, wInF);
		mpfr_set(x, m_minInputX, MPFR_RNDU);
		
		TestCase *tc=new TestCase(this);
		tc->addFPInput("X", 0.0);
		emulate(tc);
		tcl->add(tc);
		
		tc=new TestCase(this);
		FPNumber fpx(wInE,wInF,x);
		tc->addFPInput("X", &fpx);
		emulate(tc);
		tcl->add(tc);
		
		tc=new TestCase(this);
		tc->addFPInput("X", -0.125);
		emulate(tc);
		tcl->add(tc);
		
		tc=new TestCase(this);
		tc->addFPInput("X", -0.25);
		emulate(tc);
		tcl->add(tc);
		
		tc=new TestCase(this);
		tc->addFPInput("X", -0.5);
		emulate(tc);
		tcl->add(tc);
		
		tc=new TestCase(this);
		tc->addFPInput("X", -1.0);
		emulate(tc);
		tcl->add(tc);
		
		tc=new TestCase(this);
		tc->addFPInput("X", -2.0);
		emulate(tc);
		tcl->add(tc);
		
		tc=new TestCase(this);
		tc->addFPInput("X", -4.0);
		emulate(tc);
		tcl->add(tc);
		
		while(mpfr_cmp_d(x, 0)<0){
			tc=new TestCase(this);
			tc->addFPInput("X", mpfr_get_d(x,MPFR_RNDN));
			emulate(tc);
			tcl->add(tc);
			
			tc=new TestCase(this);
			tc->addFPInput("X", 0.0);
			emulate(tc);
			tcl->add(tc);
			
			tc=new TestCase(this);
			tc->addFPInput("X", 0.0);
			emulate(tc);
			tcl->add(tc);
			
			mpfr_add_d(x, x, 0.03, MPFR_RNDN);
		}
		
		mpfr_clear(x);
	}
	
	std::vector<mpz_class> FPNormalCDF::emulate_Fx(mpz_class fixX)
	{
		TestCase *tc=new TestCase(opFx);
		tc->addInput("X", fixX);
		opFx->emulate(tc);
		return tc->getExpectedOutputValues("R");
	}
	
	void FPNormalCDF::emulate_Fx(mpfr_t r, mpfr_t x)
	{
		mpfr_t tmp;
		mpfr_init2(tmp, getToolPrecision());
		mpfr_set(tmp, x, MPFR_RNDN);
		mpfr_abs(tmp, tmp, MPFR_RNDN);
		mpfr_mul_2si(tmp, tmp, -fixXMSB, MPFR_RNDN);
		REPORT(DEBUG, "  emulate_Fx: before convert to fix, x="<<mpfr_get_d(tmp, MPFR_RNDN));
		mpfr_mul_2si(tmp, tmp, (fixXMSB-fixXLSB+1), MPFR_RNDN);
		REPORT(DEBUG, "  emulate_Fx: after convert to fix, x="<<mpfr_get_d(tmp, MPFR_RNDN));
		
		mpz_class fx;
		mpfr_get_z(fx.get_mpz_t(), tmp, MPFR_RNDN);
		
		std::vector<mpz_class> v=emulate_Fx(fx);
		
		REPORT(DEBUG, "  emulate_Fx: result before convert from fix = "<<v[0]);
		
		mpfr_set_z(tmp, v[0].get_mpz_t(), MPFR_RNDN);
		mpfr_mul_2si(tmp, tmp, FxLSB, MPFR_RNDN);
		mpfr_set(r, tmp, MPFR_RNDN);
		
		mpfr_clear(tmp);
	}


	void FPNormalCDF::emulate(TestCase * tc)
	{
		std::cerr.precision(10);
		
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");

		/* Compute correct value */
		FPNumber fpx(wInE, wInF);
		fpx = svX;
		
		mpfr_t x, w, tmp, r, vUp, vDown;
		mpfr_init2(x, 1+wInF);
		fpx.getMPFR(x);
		
		mpfr_init2(w, getToolPrecision());
		mpfr_init2(tmp, getToolPrecision());
		NCD(w, x);
		
		REPORT(DEBUG, "NCD("<<mpfr_get_d(x,MPFR_RNDN)<<") = "<<mpfr_get_d(w, MPFR_RNDN));
		
		mpfr_init2(r, 1+wOutF); 
		
		mpfr_init2(vUp, 1+wOutF); 
		mpfr_set(vUp, w, MPFR_RNDU);
		REPORT(DEBUG, "  R(up="<<mpfr_get_d(vUp,MPFR_RNDN));
		FPNumber  fpru(wOutE, wOutF, vUp);
		fpru.getMPFR(vUp);
		mpz_class svru = fpru.getSignalValue();
		REPORT(DEBUG, "  R(up="<<mpfr_get_d(vUp,MPFR_RNDN)<<",  R(up,raw) = "<<svru);
		tc->addExpectedOutput("R", svru);
		
		/*if(mpfr_cmp(vUp, w) < 0){
			std::cerr<<"Ordering not maintained.\n";
			exit(1);
		}*/

		mpfr_init2(vDown, 1+wOutF);
		mpfr_set(vDown, w, MPFR_RNDD);
		REPORT(DEBUG, "  R(down="<<mpfr_get_d(vDown,MPFR_RNDN));
		mpfr_fprintf(stderr, " in = %Rb\n", vDown);
		FPNumber  fprd(wOutE, wOutF, vDown);
		mpfr_fprintf(stderr, "        %Zx\n", fprd.getMantissaSignalValue().get_mpz_t());
		fprd.getMPFR(vDown);
		mpfr_fprintf(stderr, "out = %Rb\n", vDown);
		mpz_class svrd = fprd.getSignalValue();
		tc->addExpectedOutput("R", svrd);
		REPORT(DEBUG, "  R(down="<<mpfr_get_d(vDown,MPFR_RNDN)<<",  R(down,raw) = "<<svrd);
		
		/*if(mpfr_cmp(vDown, w) > 0){
			std::cerr<<"Ordering not maintained.\n";
			exit(1);
		}*/
		
		if(m_debugOutputs){
			mpfr_prec_round(w, 1+wSqrXF, MPFR_RNDN);
			
			mpfr_t BxUp, BxDown, FxUp, FxDown, Res;
			mpfr_inits2(1+wBxF, BxUp, BxDown, NULL);
			mpfr_inits2(1+wFxF, FxUp, FxDown, NULL);
			mpfr_init2(Res, 1+wOutF);
			
			mpfr_sqr(w, x, MPFR_RNDN);
			FPNumber fpx(wInE, wSqrXF, w);
			tc->addExpectedOutput("oDebug_sqrX", fpx.getSignalValue());
			//tc->addExpectedOutput("oDebug_sqrX_expnt", fpx.getExponentSignalValue());
			//tc->addExpectedOutput("oDebug_sqrX_frac", fpx.getMantissaSignalValue());
			
			mpfr_div_si(w, w, -2, MPFR_RNDN);
			fpx=FPNumber(wInE, wSqrXF, w);
			tc->addExpectedOutput("oDebug_negHalfSqrX", fpx.getSignalValue());
			//tc->addExpectedOutput("oDebug_negHalfSqrX_expnt", fpx.getExponentSignalValue());
			//tc->addExpectedOutput("oDebug_negHalfSqrX_frac", fpx.getMantissaSignalValue());
			
			mpfr_prec_round(tmp, 1+wBxF, MPFR_RNDN);
			
			mpfr_exp(BxDown, w, MPFR_RNDD);
			fpx=FPNumber(wBxE, wBxF, BxDown);
			REPORT(FULL, "    oDebug_Bx up = "<<BxDown<<", raw="<<fpx.getSignalValue());
			tc->addExpectedOutput("oDebug_Bx", fpx.getSignalValue());
			//tc->addExpectedOutput("oDebug_Bx_expnt", fpx.getExponentSignalValue());
			//tc->addExpectedOutput("oDebug_Bx_frac", fpx.getMantissaSignalValue());
			
			mpfr_exp(BxUp, w, MPFR_RNDU);
			fpx=FPNumber(wBxE, wBxF, BxUp);
			REPORT(FULL, "    oDebug_Bx down = "<<BxUp<<", raw="<<fpx.getSignalValue());
			tc->addExpectedOutput("oDebug_Bx", fpx.getSignalValue());
			//tc->addExpectedOutput("oDebug_Bx_expnt", fpx.getExponentSignalValue());
			//tc->addExpectedOutput("oDebug_Bx_frac", fpx.getMantissaSignalValue());
			
			//////////////////////////////////////////////
			// Now do FX branch
			
			mpfr_prec_round(w, getToolPrecision(), MPFR_RNDN);
			mpfr_set(w, x, MPFR_RNDN);
			mpfr_abs(w,w,MPFR_RNDN);
			mpfr_mul_2si(w, w, -fixXLSB+1, MPFR_RNDN);
			mpz_class fixX;
			mpfr_get_z(fixX.get_mpz_t(), w, MPFR_RNDN);
			// Clamp to valid range
			fixX=std::max(mpz_class(0), std::min(mpz_class((mpz_class(1)<<(fixXMSB-fixXLSB+1))-1), fixX));
			REPORT(DEBUG, "  oDebug_fixX="<<fixX);
			tc->addExpectedOutput("oDebug_fixX", fixX);
			
			std::vector<mpz_class> FxOutputs=emulate_Fx(fixX);
			
			mpfr_div_2si(w, w, -fixXLSB+1, MPFR_RNDN);
			mpfr_prec_round(tmp, getToolPrecision(), MPFR_RNDN);
			
			mpfr_neg(w, w, MPFR_RNDN);
			F(tmp, w, MPFR_RNDN);
			
			REPORT(DEBUG, "  F("<<mpfr_get_d(w,MPFR_RNDN)<<") = "<<mpfr_get_d(tmp,MPFR_RNDN));
			{
				mpfr_t t2;
				mpfr_init2(t2,getToolPrecision());
				emulate_Fx(t2, w);
				double ratio=mpfr_get_d(tmp, MPFR_RNDN) / mpfr_get_d(t2,MPFR_RNDN);
				REPORT(DEBUG, "   emulate_fx = "<<mpfr_get_d(t2, MPFR_RNDN)<<"  ratio="<<ratio);
				
				mpfr_set(t2, w, MPFR_RNDN);
				mpfr_neg(t2, t2, MPFR_RNDN);	
				REPORT(DEBUG, "  fixXMSB="<<fixXMSB<<", multiplying by 2^"<<-fixXMSB);
				mpfr_mul_2si(t2, t2, -fixXMSB, MPFR_RNDN);
				REPORT(DEBUG, "   mapped to [0,1) = "<<mpfr_get_d(t2, MPFR_RNDN));
				funcFx->eval(t2, t2);
				ratio=mpfr_get_d(tmp, MPFR_RNDN) / mpfr_get_d(t2,MPFR_RNDN);
				REPORT(DEBUG, "   evalute_FuncStr = "<<mpfr_get_d(t2, MPFR_RNDN)<<"  ratio="<<ratio);
				mpfr_clear(t2);
				
				
			}
			
			// Convert to fixed-point
			mpfr_mul_2si(tmp, tmp, -FxLSB, MPFR_RNDN);
			
			// We have exact value in tmp, now do round up and down to fixed-point values
			mpfr_ceil(FxUp, tmp);
			mpfr_get_z(fixX.get_mpz_t(), FxUp, MPFR_RNDN);
			//tc->addExpectedOutput("oDebug_fixFx_trunc", fixX);
			REPORT(DEBUG, "   roundDownFix="<<fixX<<" vs "<<FxOutputs.front());
			mpfr_mul_2si(w, w, FxLSB, MPFR_RNDN);
			fpx=FPNumber(wFxE,wFxF, w);
			//tc->addExpectedOutput("oDebug_Fx", fpx.getSignalValue());
			mpfr_mul_2si(FxUp, FxUp, FxLSB, MPFR_RNDN);
			
			mpfr_floor(FxDown, tmp);
			mpfr_get_z(fixX.get_mpz_t(), FxDown, MPFR_RNDN);
			//tc->addExpectedOutput("oDebug_fixFx_trunc", fixX);
			REPORT(DEBUG, "   roundUpFix="<<fixX<<" vs "<<FxOutputs.back());
			mpfr_mul_2si(w, w, FxLSB, MPFR_RNDN);
			fpx=FPNumber(wFxE,wFxF, w);
			//tc->addExpectedOutput("oDebug_Fx", fpx.getSignalValue());
			mpfr_mul_2si(FxDown, FxDown, FxLSB, MPFR_RNDN);
			
			for(unsigned iF=0;iF<2;iF++){
				for(unsigned iB=0;iB<2;iB++){
					mpfr_mul(Res, iF?FxUp:FxDown, iB?BxUp:BxDown, MPFR_RNDN);
					FPNumber grr(wOutE, wOutF, Res);
					grr.getMPFR(Res);
					REPORT(DEBUG, "  possible output="<<mpfr_get_d(Res, MPFR_RNDN));
					if(mpfr_cmp(Res, vUp)!=0 && mpfr_cmp(Res, vDown)!=0){
						std::cerr<<"Doesn't match either output.\n";
						//exit(1);
					}
				}
			}
		
			mpfr_clears(BxUp, BxDown, FxUp, FxDown,NULL);	
		}

		mpfr_clears(x, r, w, tmp, vUp, vDown, NULL);
	}

}
