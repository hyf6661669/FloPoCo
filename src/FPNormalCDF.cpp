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
#include "FPAdderSinglePath.hpp"

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
		mpfr_init2(minInputX, 1+wInF);
		minInput.getMPFR(minInputX);
		REPORT(DETAILED, "minInput from input range="<<mpfr_get_d(minInputX, MPFR_RNDN));
		
		mpfr_t minUserX;
		mpfr_init2(minUserX, 1+wInF);
		if(minInputValue>=0){
			mpfr_set_d(minUserX, -8, MPFR_RNDD);
		}else{
			mpfr_set_d(minUserX, minInputValue, MPFR_RNDU);
		}
		mpfr_nextabove(minUserX);	// It is much more efficient for function eval if we are just above whatever the user wants, e.g. if it is 8
		
		if(mpfr_cmp(minInputX, minUserX)<0){
			mpfr_set(minInputX, minUserX, MPFR_RNDN);
			REPORT(DETAILED, "minInput updated from user="<<mpfr_get_d(minInputX, MPFR_RNDN));	
		}
		
		mpfr_clear(minUserX);
		
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

		mpfr_init2(m_minInputX, 1+wInF);
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

	FPNormalCDF::FPNormalCDF(Target* target, int wInE, int wInF, int wOutE, int wOutF, double minInputValue, bool debugOutputs) :
		Operator(target),
		wInE(wInE), wInF(wInF),
		wOutE(wOutE), wOutF(wOutF),
		minInputValue(minInputValue),
		funcFx(NULL),
		m_debugOutputs(debugOutputs)
	{
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
		
			F(x) = exp( -x^2*(1+epsSqr)/2) *(1+expExp)
			NCDa(x) = (1+epsMul)*F(x) * (B(x)+epsBx)
			
			There are pretty much two parts: before the exp, and after the exp. So it
			reduces to:
			
			F(x) = exp( -x^2*(1+epsPre)/2)*(1+epsPost)
			NCDa(x) = (1+epsPost)*F(x) * (B(x)+epsPost)
			
			or:
			
			F(x) = exp( -x^2*(1+epsPre)/2)
			NCDa(x) = (1+epsPost)^3 * F(x) * B(x)
			
			We only care about x<=0, as the other side we get by reflection. The output side is easy,
			it is just (1+epsPost)^3 < 4*epsPost. The error on the pre exp part is not too bad either,
			we just need, if we look at abs((F(x)-F'(x))/F(x)), then the largest error is simply for
			negative x.
			
			So for taking epsPre=2^-fPre and epsPost=2^-fPost, we can check that it satisfies a
			given output target by ensuring:
			
			2^(-wOutF-1) >= 4*2^-fPost + abs(  (exp(-minX^2/2 * (1-2^-fPre)) - exp(-minX^2/2) )/exp(-minX^2/2)) );
			
			So the strategy is to start at fPost= wOutF+3, and try to walk them up. Due to
			restrictions on exp, we need to have fPre>=fPost.
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
		vhdl<<tab<<declare("sqrX_sign")<<" <= sqrX("<<wSqrX-3<<");\n";
		vhdl<<tab<<declare("sqrX_expnt", wInE)<<" <= sqrX"<<range(wInE+wSqrXF-1,wSqrXF)<<";\n";
		vhdl<<tab<<declare("sqrX_frac", wSqrXF)<<" <= sqrX"<<range(wSqrXF-1,0)<<";\n";
		//vhdl<<tab<<declare("negHalfSqrX",wSqrX)<<" <= sqrX when not sqrX_flags=\"01\" else sqrX_flags & '1' & (sqrX_expnt-1) & sqrX_frac;\n";
		vhdl<<tab<<declare("negHalfSqrX",wSqrX)<<" <= sqrX_flags & ('1' or sqrX_sign) & (sqrX_expnt-1) & sqrX_frac when (sqrX_expnt>0)\n"; // No underflow
		vhdl<<tab<<tab<<" else sqrX_flags & '1' & "<<zg(wInE)<<" & "<<zg(wSqrXF)<<";\n";
		syncCycleFromSignal("negHalfSqrX");
		
		// Bx=exp(-abs(x)^2/2)
		wBxF=wSqrXF;
		wBxE=wInE;
		int wB=2+1+wBxE+wBxF;
		//target->setNotPipelined();
		FPExp *opExp=new FPExp(target, wInE, wSqrXF, 0, 0);
		//target->setPipelined();
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
		vhdl<<tab<<declare("reflect_at_end")<<" <= not X("<<wInE+wInF<<");\n";
		vhdl<<tab<<declare("absX",2+1+wInE+wInF)<<" <= X"<<range(2+1+wInE+wInF-1,1+wInE+wInF)<<"&('0' and X("<<(wInE+wInF)<<"))&X"<<range(wInE+wInF-1,0)<<";\n";
		
		syncCycleFromSignal("absX");
		
		// Convert it to fixed-point, and in parallel work out whether we are underlflowing
		fixXLSB=-(wInF+8);
		fixXMSB=CalcFxInputMSB()+1;
		REPORT(INFO, "fixX: MSB=2^"<<fixXMSB<<", fixXLSB=2^"<<fixXLSB);
		FP2FixV2 *opFP2Fix=new FP2FixV2(target, /*LSB0*/ fixXLSB-1, /*MSBO*/ fixXMSB-1, /*Signed*/ false, /*wER*/ wInE, wInF, /*trunc_p*/ false);
		oplist.push_back(opFP2Fix);
		inPortMap(opFP2Fix, "I", "absX");
		outPortMap(opFP2Fix, "O", "fixX");
		vhdl<<tab<<instance(opFP2Fix, "fp2fix")<<"\n";

		vhdl<<tab<<declare("have_underflow")<<" <= '1' when (absX >= \"" << fp2bin(m_minInputX, wInE, wInF) <<"\") else '0';\n";
		
		syncCycleFromSignal("fixX");
		syncCycleFromSignal("have_underflow");
		
		// Do the actual evaluation of F(x)
		FxLSB=CalcFxOutputLSB(wOutF+2);
		FxMSB=-1;	// The maximum is always 0.5
		int degree=4;	// TODO
		if(wOutF>=40){
			degree=5;
		}
		if(wOutF>=53){
			degree=7;
		}
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
		outPortMap(opFix2FP, "O", "Fx_preUnderflow");
		vhdl<<tab<<instance(opFix2FP, "fix2fp")<<"\n";
		syncCycleFromSignal("Fx_preUnderflow");
		
		vhdl<<tab<<declare("Fx",3+wFxE+wFxF)<<" <= "<< zg(3+wFxE+wFxF)<< " when have_underflow='1' else Fx_preUnderflow;\n";
		
		
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
		outPortMap(opMult, "R", "neg_result");
		vhdl<<tab<<instance(opMult, "finalMult")<<"\n";
		syncCycleFromSignal("neg_result");
		
		vhdl<<tab<<declare("minus_neg_result",3+wOutE+wOutF)<<" <= ";
		vhdl<<"neg_result"<<range(2+wOutE+wOutF,1+wOutE+wOutF)<<" & '1' & neg_result"<<range(wOutE+wOutF-1,0)<<";\n";
		
		// This is horribly wasteful, just to do 1-x
		int wOneE=wOutE, wOneF=wOutF;	// Current adders don't allow custom widths.
		Operator *finalSub=new FPAdderSinglePath(target, /*wEX*/wOneE, /*wFX*/ wOneF, /*wEY*/ wOutE, /*wFY*/wOutF,/*wER*/ wOutE, /*wFR*/ wOutF);
		oplist.push_back(finalSub);
		{	// Pure laziness to get a constant one
			mpfr_t vOne;
			mpfr_init2(vOne, 1+wOneF);
			mpfr_set_d(vOne, 1.0, MPFR_RNDN);
			inPortMapCst(finalSub, "X", std::string("\"")+fp2bin(vOne, wOneE, wOneF)+"\"");	// Constant one
			mpfr_clear(vOne);
		}
		inPortMap(finalSub, "Y", "minus_neg_result");
		outPortMap(finalSub, "R", "pos_result");
		vhdl<<instance(finalSub, "finalReflect");
		syncCycleFromSignal("pos_result");
		
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
		
		vhdl<<tab<<"R <= pos_result when reflect_at_end='1' else neg_result;\n";
		
		REPORT(DETAILED, "Done.");
	}
  
  FPNormalCDF::~FPNormalCDF() {
	  mpfr_clear(m_minInputX);
	  if(funcFx)
		delete funcFx;
  }
  
	void FPNormalCDF::buildStandardTestCases(TestCaseList* tcl)
	{
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
		
		for(int i=FPNumber::plusInfty;i<FPNumber::smallestNegative+1;i++){
			tc=new TestCase(this);
			fpx=FPNumber(wInE,wInF,static_cast<FPNumber::SpecialValue>(i));
			tc->addFPInput("X", &fpx);
			emulate(tc);
			tcl->add(tc);
		}
		
		// Go either side of the underflow boundary by a few ULPs
		// First go below
		mpfr_set(x, m_minInputX, MPFR_RNDU);
		for(unsigned i=0;i<16;i++){
			tc=new TestCase(this);
			tc->addFPInput("X", mpfr_get_d(x,MPFR_RNDN));
			emulate(tc);
			tcl->add(tc);
			
			mpfr_neg(x,x,MPFR_RNDN);
			tc=new TestCase(this);
			tc->addFPInput("X", mpfr_get_d(x,MPFR_RNDN));
			emulate(tc);
			tcl->add(tc);
			mpfr_neg(x,x,MPFR_RNDN);
			
			mpfr_nextbelow(x);
		}
		// Then above
		mpfr_set(x, m_minInputX, MPFR_RNDU);
		for(unsigned i=0;i<8;i++){
			tc=new TestCase(this);
			tc->addFPInput("X", mpfr_get_d(x,MPFR_RNDN));
			emulate(tc);
			tcl->add(tc);
			
			mpfr_neg(x,x,MPFR_RNDN);
			tc=new TestCase(this);
			tc->addFPInput("X", mpfr_get_d(x,MPFR_RNDN));
			emulate(tc);
			tcl->add(tc);
			mpfr_neg(x,x,MPFR_RNDN);
			
			mpfr_nextabove(x);
		}
		
		
		// Now walk over values close to zero, as that is another danger zone (thought
		// for this specific method, probably not)
		mpfr_set_d(x, 0.0, MPFR_RNDU);
		for(unsigned i=0;i<8;i++){
			tc=new TestCase(this);
			tc->addFPInput("X", mpfr_get_d(x,MPFR_RNDN));
			emulate(tc);
			tcl->add(tc);
			
			mpfr_neg(x,x,MPFR_RNDN);
			tc=new TestCase(this);
			tc->addFPInput("X", mpfr_get_d(x,MPFR_RNDN));
			emulate(tc);
			tcl->add(tc);
			mpfr_neg(x,x,MPFR_RNDN);
			
			mpfr_nextbelow(x);
		}
		
		// Now just choose values equally spread across the range.
		mpfr_set(x, m_minInputX, MPFR_RNDU);
		mpfr_mul_d(x, x, 1.1, MPFR_RNDN);	// Go slightly outside standard range
		while(mpfr_cmp_d(x, 0)<0){
			tc=new TestCase(this);
			tc->addFPInput("X", mpfr_get_d(x,MPFR_RNDN));
			emulate(tc);
			tcl->add(tc);
			
			mpfr_neg(x,x,MPFR_RNDN);
			tc=new TestCase(this);
			tc->addFPInput("X", mpfr_get_d(x,MPFR_RNDN));
			emulate(tc);
			tcl->add(tc);
			mpfr_neg(x,x,MPFR_RNDN);
			
			mpfr_add_d(x, x, 0.3, MPFR_RNDN);
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
		
		mpfr_init2(w, 10+wOutF);	// NCD is well behaved, but slow
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
		

		mpfr_init2(vDown, 1+wOutF);
		mpfr_set(vDown, w, MPFR_RNDD);
		REPORT(DEBUG, "  R(down="<<mpfr_get_d(vDown,MPFR_RNDN));
		//mpfr_fprintf(stderr, " in = %Rb\n", vDown);
		FPNumber  fprd(wOutE, wOutF, vDown);
		//mpfr_fprintf(stderr, "        %Zx\n", fprd.getMantissaSignalValue().get_mpz_t());
		fprd.getMPFR(vDown);
		//mpfr_fprintf(stderr, "out = %Rb\n", vDown);
		mpz_class svrd = fprd.getSignalValue();
		tc->addExpectedOutput("R", svrd);
		REPORT(DEBUG, "  R(down="<<mpfr_get_d(vDown,MPFR_RNDN)<<",  R(down,raw) = "<<svrd);

		
		if(m_debugOutputs){
			mpfr_t BxUp, BxDown, FxUp, FxDown, Res;
			
			mpfr_inits2(1+wBxF, BxUp, BxDown, NULL);
			mpfr_inits2(1+wFxF, FxUp, FxDown, NULL);
			mpfr_init2(Res, 1+wOutF);
			
			mpfr_prec_round(w, 1+wSqrXF, MPFR_RNDN);
			
			mpfr_set(tmp, x, MPFR_RNDN);
			mpfr_abs(tmp, tmp, MPFR_RNDN);
			mpfr_neg(tmp, tmp, MPFR_RNDN);
			
			if(mpfr_cmp(tmp, m_minInputX) > 0){
				
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
				
				/*
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
				*/
			}
		
			mpfr_clears(BxUp, BxDown, FxUp, FxDown,NULL);	
		}

		mpfr_clears(x, r, w, tmp, vUp, vDown, NULL);
	}
	
	TestCase* FPNormalCDF::buildRandomTestCase(int i)
	{
		if(getLargeRandom(3)==0){
			return Operator::buildRandomTestCase(i);
		}else{
			mpfr_t input;
			mpfr_init2(input, getToolPrecision());
			mpfr_urandomb(input, FloPoCoRandomState::m_state);
			mpfr_mul(input, input, m_minInputX, MPFR_RNDN);
			mpfr_mul_d(input,input,1.05,MPFR_RNDN);
			
			if(getLargeRandom(1)==0){
				mpfr_neg(input, input, MPFR_RNDN);
			}
				
			
			TestCase *tc=new TestCase(this);
			FPNumber fpx(wInE,wInF,input);
			tc->addFPInput("X", &fpx);
			emulate(tc);
			
			mpfr_clear(input);
			
			return tc;
		}
	}

}
