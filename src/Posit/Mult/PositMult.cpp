/*
  A Posit Multiplier for FloPoCo

  Authors: Raul Murillo, Alberto A. del Barrio, Guillermo Botella

  This file is part of the research article
  "Template-Based Posit Multiplication for Training and Inferring in Neural Networks"
  (https://arxiv.org/abs/1907.04091)

  Initial software.
  Copyright © Complutense University of Madrid, Spain, 2019.

  All rights reserved.

*/
#include <sstream>
#include <iostream>
#include "PositMult.hpp"
#include <Operator.hpp>
#include <utils.hpp>

#include "PositDecoder.hpp"
#include "ShiftersEtc/Shifters.hpp"
#include "IntMult/IntMultiplier.hpp"
#include "TestBenches/PositNumber.hpp"

using namespace std;
namespace flopoco {


	PositMult::PositMult(Target* target, Operator* parentOp, int width_, int wES_):	Operator(parentOp, target), width(width_), wES(wES_) {
		srcFileName="PositMult";

		ostringstream name;
		name << "PositMult_" << width << "_" << wES ;
		setNameWithFreqAndUID(name.str());
		setCopyrightString("Raul Murillo, Alberto A. del Barrio, Guillermo Botella, 2019");

		if (width < 3) {
			throw std::string("PositMult Constructor : width is too small, should be greater than two");
		}
		int freeWidth = width - 3;
		if (wES >= freeWidth) {
			//Avoid posits without even one bit of precision
			throw std::string("PositMult Constructor : invalid value of wES");
		}

		// SET UP THE IO SIGNALS

		// declaring inputs
		int RegSize = intlog2(width);
		int FracSize = width-wES-2;
		addInput  ( "X", width);
		addInput  ( "Y", width);
		// declaring output
		addOutput ( "R", width);

		addFullComment("Start of vhdl generation");

		REPORT(INFO,"Declaration of PositMult \n");
		REPORT(DETAILED, "this operator has received two parameters " << width << " and " << wES);
  		REPORT(DEBUG,"debug of PositMult");

	//=========================================================================|
		addFullComment("Data Extraction");
	// ========================================================================|

		ostringstream paramX, inmapX, outmapX;
		paramX << "width=" << width;
		paramX << " wES=" << wES;
		
		inmapX << "Input=>X";
		
		outmapX << "Sign=>sign_X,Reg=>reg_X,Exp=>exp_X,Frac=>frac_X,z=>z_X,inf=>inf_X";
		
		newInstance("PositDecoder", "X_decoder", paramX.str(), inmapX.str(), outmapX.str());

		ostringstream paramY, inmapY, outmapY;
		paramY << "width=" << width;
		paramY << " wES=" << wES;
		
		inmapY << "Input=>Y";
		
		outmapY << "Sign=>sign_Y,Reg=>reg_Y,Exp=>exp_Y,Frac=>frac_Y,z=>z_Y,inf=>inf_Y";
		
		newInstance("PositDecoder", "Y_decoder", paramY.str(), inmapY.str(), outmapY.str());

		addComment("Gather scale factors");
		vhdl << declare(.0, "sf_X", RegSize+wES) << " <= reg_X";
		if (wES>0){
			vhdl << " & exp_X";
		} vhdl << ";" << endl;
		vhdl << declare(.0, "sf_Y", RegSize+wES) << " <= reg_Y";
		if (wES>0){
			vhdl << " & exp_Y";
		} vhdl << ";" << endl;
		
	//=========================================================================|
		addFullComment("Sign and Special Cases Computation");
	// ========================================================================|
	
		vhdl << declare(target->logicDelay(1), "sign", 1, false) << " <= sign_X XOR sign_Y;" << endl;
		vhdl << declare(target->logicDelay(1), "z", 1, false) << " <= z_X OR z_Y;" << endl;
		vhdl << declare(target->logicDelay(1), "inf", 1, false) << " <= inf_X OR inf_Y;" << endl;

	//=========================================================================|
		addFullComment("Multiply the fractions, add the exponent values");
	// ========================================================================|

		ostringstream param_mult, inmap_mult, outmap_mult;
		param_mult << "wX=" << FracSize;
		param_mult << " wY=" << FracSize;
		//param_mult << " wOut=" << 0;
		param_mult << " dspThreshold=" << 0.0;

		inmap_mult << "X=>frac_X,Y=>frac_Y";

		outmap_mult << "R=>frac_mult";

		newInstance("IntMultiplier", "mantissa_multiplier", param_mult.str(), inmap_mult.str(), outmap_mult.str());

		int mult_size = 2 * FracSize; // This must be equal to getSignalByName("frac_mult")->width();

		addComment("Adjust for overflow");
		vhdl << declare(.0, "ovf_m", 1, false) << " <= frac_mult(frac_mult'high);" << endl;

		vhdl << "with ovf_m select " << declare(target->logicDelay(1), "normFrac", mult_size+1) << "<= " << endl <<
      		tab << "frac_mult & '0' when '0'," << endl <<
			tab << "'0' & frac_mult when '1'," << endl <<
			tab << "\"" << string(mult_size+1, '-') << "\" when others;" << endl;
			// Equivalent to shift right ovf_m bits

		vhdl << declare(target->adderDelay(RegSize+wES+1), "sf_mult", RegSize+wES+1) << " <= (sf_X(sf_X'high) & sf_X) + (sf_Y(sf_Y'high) & sf_Y) + ovf_m;" << endl;
		vhdl << declare(.0, "sf_sign") << " <= sf_mult(sf_mult'high);" << endl;

	//=========================================================================|
		addFullComment("Compute Regime and Exponent value");
	// ========================================================================|

		vhdl << declare(target->adderDelay(mult_size), "nzero", 1, false) << " <= '0' when frac_mult = \"" << string(mult_size, '0') << "\" else '1';"<<endl;		

		addComment("Unpack scaling factors");
		if (wES > 0){
			vhdl << declare(.0, "ExpBits", wES) << " <= sf_mult" << range(wES - 1, 0) << ";" << endl;
		}
		vhdl << declare(.0, "RegimeAns_tmp", RegSize) << " <= sf_mult" << range(RegSize + wES - 1, wES) << ";" << endl;

		addComment("Get Regime's absolute value");
		vhdl << "with sf_sign select " << declare(target->logicDelay() + target->adderDelay(RegSize), "RegimeAns", RegSize) << "<=" << endl <<
			tab << "(NOT RegimeAns_tmp)+1 when '1'," << endl <<
			tab << "RegimeAns_tmp when '0'," << endl <<
			tab << "\"" << string(RegSize, '-') << "\" when others;" << endl;

		addComment("Check for Regime overflow");
		vhdl << declare(.0, "ovf_reg", 1, false) << " <= RegimeAns(RegimeAns'high);" << endl;
		vhdl << "with ovf_reg select " << declare(target->logicDelay(1), "FinalRegime", RegSize) << "<=" << endl <<
			tab << "'0' & \"" << string(RegSize-1, '1') << "\" when '1'," << endl <<
      		tab << "RegimeAns when '0'," << endl <<
      		tab << "\"" << string(RegSize, '-') << "\" when others;" << endl;

		vhdl << declare(target->adderDelay(RegSize), "ovf_regF", 1, false) << " <= '1' when FinalRegime = ('0' & \"" << string(RegSize-1, '1') << "\") else '0';" << endl;		
		vhdl << declare(target->logicDelay(1), "ovf_anyreg", 1, false) << " <= ovf_reg OR ovf_regF;" << endl;		
		if (wES > 0){
			vhdl << declare(target->logicDelay(2), "ovf_exp", 1, false) << " <= ovf_anyreg OR (NOT nzero);" << endl;
			vhdl << "with ovf_exp select " << declare(target->logicDelay(1), "FinalExp", wES) << "<=" << endl <<
				tab << "\"" << string(wES, '0') << "\" when '1'," << endl <<
				tab << "ExpBits when '0'," << endl <<
      			tab << "\"" << string(wES, '-') << "\" when others;" << endl;
		}
		
	//=========================================================================|
		addFullComment("Packing Stage 1");
	// ========================================================================|

		vhdl << declare(.0, "tmp1", 2+wES+mult_size-1) << " <= nzero & '0' ";
		if (wES > 0){
			vhdl << "& FinalExp ";
		}
		vhdl << "& normFrac" << range(mult_size-2, 0) << ";" << endl;
		vhdl << declare(.0, "tmp2", 2+wES+mult_size-1) << " <= '0' & nzero ";
		if (wES > 0){
			vhdl << "& FinalExp ";
		}
		vhdl << "& normFrac" << range(mult_size-2, 0) << ";" << endl;

		vhdl << "with ovf_regF select " << declare(target->adderDelay(RegSize), "shift_neg", RegSize) << "<=" << endl <<
			tab << "FinalRegime - 2 when '1'," << endl <<
			tab << "FinalRegime - 1 when '0'," << endl <<
			tab << "\"" << string(RegSize, '-') << "\" when others;" << endl;
		vhdl << "with ovf_regF select " << declare(target->adderDelay(RegSize), "shift_pos", RegSize) << "<=" << endl <<
			tab << "FinalRegime - 1 when '1'," << endl <<
			tab << "FinalRegime when '0'," << endl <<
			tab << "\"" << string(RegSize, '-') << "\" when others;" << endl;

		vhdl << "with sf_sign select " << declare(target->logicDelay(1), "input_shifter", 2+wES+mult_size-1) << "<=" << endl <<
			tab << "tmp2 when '1'," << endl <<
			tab << "tmp1 when '0'," << endl <<
			tab << "\"" << string(2+wES+mult_size-1, '-') << "\" when others;" << endl;
		vhdl << "with sf_sign select " << declare(target->logicDelay(1), "shift_offset", RegSize) << "<=" << endl <<
			tab << "shift_neg when '1'," << endl <<
			tab << "shift_pos when '0'," << endl <<
			tab << "\"" << string(RegSize, '-') << "\" when others;" << endl;
		vhdl << declare(.0, "pad", 1, false) << "<= input_shifter(input_shifter'high);" << endl;

		ostringstream param_shift, inmap_shift, outmap_shift;
		param_shift << "wIn=" << 2 + wES + mult_size - 1;
		param_shift << " maxShift=" << width;
		param_shift << " wOut=" << width + 1;	// width-1 + (G, R bits) PROBLEM: if(computeSticky) wOut =  wIn;
		param_shift << " dir=" << Shifter::Right;
		param_shift << " computeSticky=true";
		param_shift << " inputPadBit=true";
		
		inmap_shift << "X=>input_shifter,S=>shift_offset,padBit=>pad";
		
		outmap_shift << "R=>shifted_frac,Sticky=>S_bit_tmp";
		
		newInstance("Shifter", "right_signed_shifter", param_shift.str(), inmap_shift.str(), outmap_shift.str());
		int shift_size = getSignalByName("shifted_frac")->width();

		vhdl << declare(.0, "tmp_ans", width-1) << " <= shifted_frac" << range(shift_size-1, shift_size-(width-1)) << ";" << endl;

	//=========================================================================|
		addFullComment("Packing Stage 2 - Unbiased Rounding");
	// ========================================================================|
		// Rounding implementation using L,G,R,S bits
		vhdl << declare(.0, "LSB", 1, false) << " <= shifted_frac" << of(shift_size-(width-1)) << ";" << endl;
		vhdl << declare(.0, "G_bit", 1, false) << " <= shifted_frac" << of(shift_size-(width-1)-1) << ";" << endl;
		vhdl << declare(.0, "R_bit", 1, false) << " <= shifted_frac" << of(shift_size-(width-1)-2) << ";" << endl;
		vhdl << declare(target->adderDelay(shift_size-(width-1)-2), "S_rem", 1, false) << " <= '0' when shifted_frac" << range(shift_size-(width-1)-3, 0) << " = \"" << string(shift_size-(width-1)-2, '0') << "\" else '1';" << endl;
		vhdl << declare(target->logicDelay(1),"S_bit", 1, false) << " <= S_bit_tmp OR S_rem;" << endl;

		vhdl << "with ovf_anyreg select " << declare(target->logicDelay(3), "round", 1, false) << "<=" << endl <<
			tab << "'0' when '1'," << endl <<
			tab << "G_bit AND (LSB OR R_bit OR S_bit) when '0'," << endl <<
			tab << "'-' when others;" << endl;

		vhdl << "with sign select " << declare(target->adderDelay(width-1) + target->logicDelay() + target->adderDelay(width), "result", width) << "<=" << endl <<
			tab << "'0' & (tmp_ans + round) when '0'," << endl <<
			tab << "NOT('0' & (tmp_ans + round))+1 when '1'," << endl <<
			tab << "\"" << string(width, '-') << "\" when others;" << endl;

		vhdl << "R <= '1' & \"" << string(width-1, '0') << "\" when inf = '1' else " << endl <<
			tab << zg(width) << " when z = '1' else" << endl <<
			tab << " result;" << endl;


		addFullComment("End of vhdl generation");

	};

	
	void PositMult::emulate(TestCase * tc) {
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
		
		/* Compute correct value */
		PositNumber posx(width, wES, svX);
		PositNumber posy(width, wES, svY);
		mpfr_t x, y, r;
		mpfr_init2(x, 1000*width -2);
		mpfr_init2(y, 1000*width -2);
		mpfr_init2(r, 1000*width -2);
		posx.getMPFR(x);
		posy.getMPFR(y);
		mpfr_mul(r, x, y, GMP_RNDN);
		
		// Set outputs
		PositNumber posr(width, wES, r);
		mpz_class svR = posr.getSignalValue();
		tc->addExpectedOutput("R", svR);
		
		// clean up
		mpfr_clears(x, y, r, NULL);

	}


//	void PositMult::buildStandardTestCases(TestCaseList * tcl) {
		// please fill me with regression tests or corner case tests!
//	}


	OperatorPtr PositMult::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
		int width, wES;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width); 
    	UserInterface::parsePositiveInt(args, "wES", &wES);
    	return new PositMult(target, parentOp, width, wES);

	}

	void PositMult::registerFactory(){
		UserInterface::add("PositMult", // name
						"A posit multiplier with a single architecture.",
						"BasicPosit",
						"", //seeAlso
						"width(int): posit size in bits; \
                        				wES(int): exponent size in bits;",
						"",
						PositMult::parseArguments
						) ;
	}

}//namespace
