/*
  A Posit Adder for FloPoCo

  Authors: Raul Murillo, Alberto A. del Barrio, Guillermo Botella
  
  This file is part of the research article
  "Customized Posit Adders and Multipliers using the FloPoCo Core Generator"

  Initial software.
  Copyright © Complutense University of Madrid, Spain,
  2020.
  All rights reserved.
  
*/
#include <sstream>
#include <iostream>
#include "PositAdd2.hpp"
#include <Operator.hpp>
#include <utils.hpp>

#include "Posit/Decoder/PositDecoder.hpp"
#include "ShiftersEtc/Shifters.hpp"
#include "TestBenches/PositNumber.hpp"

using namespace std;
namespace flopoco {


	std::string toBinary(int n, int len){
	    std::string r;
	    while(n!=0) {r=(n%2==0 ? "0":"1")+r; n/=2;}
	    return "\""+string(len - r.length(), '0') + r+"\"";
	}


	PositAdd2::PositAdd2(Target* target, Operator* parentOp, int width_, int wES_):	Operator(parentOp, target), width(width_), wES(wES_) {
		srcFileName="PositAdd2";

		ostringstream name;
		name << "PositAdd2_" << width << "_" << wES ;
		setNameWithFreqAndUID(name.str());
		setCopyrightString("Raul Murillo, Alberto A. del Barrio, Guillermo Botella, 2020");

		if (width < 3) {
			throw std::string("PositAdd2 Constructor : width is too small, should be greater than two");
		}
		int freeWidth = width - 3;
		if (wES >= freeWidth) {
			//Avoid posits without even one bit of precision
			throw std::string("PositAdd2 Constructor : invalid value of wES");
		}

		// SET UP THE IO SIGNALS

		// declaring inputs
		int RegSize = intlog2(width-1)+1;
		int FracSize = width-wES-2;
		addInput  ( "X", width);
		addInput  ( "Y", width);
		// declaring output
		addOutput ( "R", width);

		addFullComment("Start of vhdl generation");

		REPORT(INFO,"Declaration of PositAdd2 \n");
		REPORT(DETAILED, "this operator has received two parameters " << width << " and " << wES);
  		REPORT(DEBUG,"debug of PositAdd2");

	//=========================================================================|
		addFullComment("Data Extraction");
	// ========================================================================|

		ostringstream paramX, inmapX, outmapX;
		paramX << "width=" << width;
		paramX << " wES=" << wES;
		
		inmapX << "Input=>X";
		
		//outmapX << "Sign=>sign_X,Reg=>reg_X,Exp=>exp_X,Frac=>a_frac_X,z=>z_X,inf=>inf_X,Abs_in=>X_abs";
		outmapX << "Sign=>sign_X,Reg=>reg_X,Exp=>exp_X,Frac=>frac_X,z=>open,inf=>inf_X,Abs_in=>X_abs";
		
		newInstance("PositDecoder", "X_decoder", paramX.str(), inmapX.str(), outmapX.str());

		ostringstream paramY, inmapY, outmapY;
		paramY << "width=" << width;
		paramY << " wES=" << wES;
		
		inmapY << "Input=>Y";
		
		//outmapY << "Sign=>sign_Y,Reg=>reg_Y,Exp=>exp_Y,Frac=>a_frac_Y,z=>z_Y,inf=>inf_Y,Abs_in=>Y_abs";
		outmapY << "Sign=>sign_Y,Reg=>reg_Y,Exp=>exp_Y,Frac=>frac_Y,z=>open,inf=>inf_Y,Abs_in=>Y_abs";
		
		newInstance("PositDecoder", "Y_decoder", paramY.str(), inmapY.str(), outmapY.str());

		/*//
		vhdl << declare(.0, "frac_X", FracSize) << " <= '1' & a_frac_X" << range(FracSize-2,0) << ";" << endl;
		vhdl << declare(.0, "frac_Y", FracSize) << " <= '1' & a_frac_Y" << range(FracSize-2,0) << ";" << endl;
		//*/
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
	
		vhdl << declare(target->logicDelay(1), "OP", 1, false) << " <= sign_X XOR sign_Y;" << endl;
		vhdl << declare(target->logicDelay(1), "inf", 1, false) << " <= inf_X OR inf_Y;" << endl;

	//=========================================================================|
		addFullComment("Compare operands and adjust values");
	// ========================================================================|
		double delayAdderIsLargeExp = target->adderDelay(width-1);
  		cerr << "AdderIsLarge delay : " << delayAdderIsLargeExp << endl; 
    	vhdl << declare(delayAdderIsLargeExp, "is_larger", 1, false) << "<= '1' when X_abs > Y_abs else '0';" << endl;
		
		vhdl << "with is_larger select " << declare(target->logicDelay(1), "larger_sign", 1, false) << " <= " << endl <<
      		tab << "sign_X when '1'," << endl <<
     		tab << "sign_Y when '0'," << endl <<
     		tab << "'-' when others;" << endl;

		vhdl << "with is_larger select " << declare(target->logicDelay(1), "larger_sf", RegSize+wES) << " <= " << endl <<
			tab << "sf_X when '1'," << endl <<
			tab << "sf_Y when '0'," << endl <<
			tab << "\"" << string(RegSize+wES, '-') << "\" when others;" << endl;
		/*
		vhdl << "with is_larger select " << declare(target->logicDelay(1), "smaller_sf", RegSize+wES) << " <= " << endl <<
			tab << "sf_Y when '1'," << endl <<
			tab << "sf_X when '0'," << endl <<
			tab << "\"" << string(RegSize+wES, '-') << "\" when others;" << endl;
		*/
		vhdl << "with is_larger select " << declare(target->logicDelay(1), "larger_frac", FracSize) << " <= " << endl <<
			tab << "frac_X when '1'," << endl <<
			tab << "frac_Y when '0'," << endl <<
			tab << "\"" << string(FracSize, '-') << "\" when others;" << endl;

		
		// Smaller fraction is negated if operators signs are not equal, in order to perform 2's comp. addition of fractions
	/*	vhdl << declare(.0, "rep_OP", FracSize+2) << " <= (others => OP);" << endl;
		vhdl << "with is_larger select " << declare(target->logicDelay() + target->adderDelay(FracSize+2), "smaller_frac", FracSize+2) << " <= " << endl <<
			tab << "(rep_OP XOR (frac_Y & \"00\")) + OP when '1'," << endl <<
			tab << "(rep_OP XOR (frac_X & \"00\")) + OP when '0'," << endl <<
			tab << "\"" << string(FracSize+2, '-') << "\" when others;" << endl;
	*/
	
		vhdl << "with is_larger select " << declare(target->logicDelay(), "smaller_frac", FracSize+2) << " <= " << endl <<
			tab << "(frac_Y & \"00\")  when '1'," << endl <<
			tab << "(frac_X & \"00\")  when '0'," << endl <<
			tab << "\"" << string(FracSize+2, '-') << "\" when others;" << endl;
				
		// Compute smaller fraction offset as the absolute value of Scaling Factors' difference, (and saturate it)
		vhdl << declare(target->adderDelay(RegSize+wES+1), "sf_diff", RegSize+wES+1) << " <= (sf_X(sf_X'high) & sf_X) - (sf_Y(sf_Y'high) & sf_Y);" << endl;
		vhdl << declare(.0, "diff_msb", RegSize+wES) << " <= (others => sf_diff(sf_diff'high));" << endl;
		vhdl << declare(target->adderDelay(RegSize+wES) + target->logicDelay(), "offset", RegSize+wES) << " <= (diff_msb XOR sf_diff" << range(RegSize+wES-1, 0) <<") + sf_diff(sf_diff'high);" << endl;
		
		int maxshiftsize = intlog2(FracSize+2);
		if ((RegSize+wES) > maxshiftsize){
			vhdl << declare(.0, "sup_offset", RegSize+wES-maxshiftsize) << " <= offset" << range(RegSize+wES-1, maxshiftsize) << ";" << endl;

			vhdl << declare(target->adderDelay(RegSize+wES-maxshiftsize), "shift_saturate", 1, false) << " <= '0' when sup_offset = \"" << string(RegSize+wES-maxshiftsize, '0') << "\" else '1';"<<endl;

			vhdl << "with shift_saturate select " << declare(target->logicDelay(1), "frac_offset", maxshiftsize) << " <=" << endl <<
				tab << "\"" << string(maxshiftsize, '1') << "\" when '1'," << endl <<
				tab << "offset" << range(maxshiftsize-1, 0) << " when '0'," << endl <<
				tab << "\"" << string(maxshiftsize, '-') << "\" when others;" << endl;
		}
		else {	// Only for wES==0. In this case (RegSize+wES) == maxshiftsize
			if((RegSize+wES) != maxshiftsize){
				cerr << "Wrong sizes: " << (RegSize+wES) << " and " << maxshiftsize << endl;
				// Unreachable point
			}
			vhdl << declare(.0, "frac_offset", maxshiftsize) << " <= offset;" << endl;
		}
		

	//=========================================================================|
		addFullComment("Align mantissas");
	// ========================================================================|
		//*****
		//vhdl << declare(.0, "offset", intlog2(FracSize+2)) << " <= (others => '0');" << endl;
		//*****
		ostringstream param1, inmap1, outmap1;
		param1 << "wIn=" << FracSize+2;
		param1 << " maxShift=" << FracSize+2;	// Adjust this. Is there a lower maximum length to shift?
		param1 << " wOut=" << FracSize+2;
		param1 << " dir=" << Shifter::Right;
		param1 << " computeSticky=true";
		//param1 << " inputPadBit=true";
    
    	//inmap1 << "X=>smaller_frac,S=>frac_offset,padBit=>OP";
		inmap1 << "X=>smaller_frac,S=>frac_offset";
    
		outmap1 << "R=>shifted_frac,Sticky=>sticky";
		
		newInstance("Shifter", "mantissa_shifter", param1.str(), inmap1.str(), outmap1.str());
	/*
		vhdl << declare(.0, "smaller_shifted_frac", FracSize) << " <= shifted_frac" << range(FracSize+1,2) << ";" << endl;
		vhdl << declare(.0, "G_aux", 1, false) << " <= shifted_frac" << of(1) << ";" << endl;
		vhdl << declare(.0, "R_aux", 1, false) << " <= shifted_frac" << of(0) << ";" << endl;
	*/
		/////
	/*	vhdl << declare(.0, "larger_op_frac", FracSize+1) << " <= '0' & larger_frac;" << endl;
		vhdl << "with z_Y select " << declare(target->logicDelay(), "smaller_op_frac", FracSize+1) << " <=" << endl <<
			tab << "OP & smaller_shifted_frac when '0'," << endl <<
			tab << "\"" << string(FracSize+1, '0') << "\" when '1'," << endl <<
			tab << "\"" << string(FracSize+1, '-') << "\" when others;" << endl;
//*/
	//=========================================================================|
		addFullComment("Add aligned mantissas");
	// ========================================================================|
		//vhdl << declare(target->adderDelay(FracSize+1), "add_frac", FracSize+1) << " <= larger_op_frac + smaller_op_frac;" << endl;
	//	vhdl << declare(target->adderDelay(FracSize+1), "add_frac", FracSize+1) << " <= ('0' & larger_frac) + (OP & smaller_shifted_frac);" << endl;
		vhdl << "with OP select " << declare(target->adderDelay(FracSize+4), "add_frac", FracSize+4) << " <= " << endl <<
			tab << "('0' & larger_frac & \"000\") + ('0' & shifted_frac & sticky) when '0'," << endl <<
			tab << "('0' & larger_frac & \"000\") - ('0' & shifted_frac & sticky) when '1'," << endl <<
			tab << "\"" << string(FracSize+4, '-') << "\" when others;" << endl;

		vhdl << declare(.0, "ovf_frac", 1, false) << " <= add_frac(add_frac'high);" << endl;
		vhdl << "with ovf_frac select " << declare(target->logicDelay(), "useful_frac", FracSize) << " <=" << endl <<
			tab << "add_frac" << range(FracSize+3, 4) << " when '1'," << endl <<
			tab << "add_frac" << range(FracSize+2, 3) << " when '0'," << endl <<
			tab << "\"" << string(FracSize, '-') << "\" when others;" << endl;
		vhdl << "with ovf_frac select " << declare(target->logicDelay(), "G_tmp", 1, false) << " <=" << endl <<
			tab << "add_frac" << of(3) << " when '1'," << endl <<
			tab << "add_frac" << of(2) << " when '0'," << endl <<
			tab << "'-' when others;" << endl;
		vhdl << "with ovf_frac select " << declare(target->logicDelay(), "R_tmp", 1, false) << " <=" << endl <<
			tab << "add_frac" << of(2) << " when '1'," << endl <<
			tab << "add_frac" << of(1) << " when '0'," << endl <<
			tab << "'-' when others;" << endl;
		vhdl << "with ovf_frac select " << declare(target->logicDelay(), "S_tmp", 1, false) << " <=" << endl <<
			tab << "add_frac" << of(1) << " OR add_frac" << of(0) << " when '1'," << endl <<
			tab << "add_frac" << of(0) << " when '0'," << endl <<
			tab << "'-' when others;" << endl;

	/*
		vhdl << declare(.0, "ovf_frac", 1, false) << " <= add_frac(add_frac'high);" << endl;
		vhdl << "with ovf_frac select " << declare(target->logicDelay(), "useful_frac", FracSize) << " <=" << endl <<
			tab << "add_frac" << range(FracSize, 1) << " when '1'," << endl <<
			tab << "add_frac" << range(FracSize-1, 0) << " when '0'," << endl <<
			tab << "\"" << string(FracSize, '-') << "\" when others;" << endl;
		//vhdl << declare(.0, "useful_frac", FracSize+2) << " <= add_frac" << range(FracSize+1, 0) << ";" << endl; // & sticky
		vhdl << "with ovf_frac select " << declare(target->logicDelay(), "G_tmp", 1, false) << " <=" << endl <<
			tab << "add_frac" << of(0) << " when '1'," << endl <<
			tab << "G_aux when '0'," << endl <<
			tab << "'-' when others;" << endl;
		vhdl << "with ovf_frac select " << declare(target->logicDelay(), "R_tmp", 1, false) << " <=" << endl <<
			tab << "G_aux when '1'," << endl <<
			tab << "R_aux when '0'," << endl <<
			tab << "'-' when others;" << endl;
		vhdl << "with ovf_frac select " << declare(target->logicDelay(), "S_tmp", 1, false) << " <=" << endl <<
			tab << "R_aux OR sticky when '1'," << endl <<
			tab << "sticky when '0'," << endl <<
			tab << "'-' when others;" << endl;

	*/

		vhdl << declare(.0, "frac_GR", FracSize+2) << " <= useful_frac & G_tmp & R_tmp;" << endl;

		addComment("Normalization of add_frac");
		ostringstream param2, inmap2, outmap2;
		int wCount = intlog2(FracSize+2);	//?? 
		param2 << "wIn=" << FracSize+2;	//??
		param2 << " wOut=" << FracSize+2;	//??
		param2 << " wCount=" << wCount;	//?? 
		param2 << " countType=" << 0;
		
		inmap2 << "I=>frac_GR";
		
		outmap2 << "Count=>lzCount,O=>normFrac";
		
		newInstance("LZOCShifterSticky", "align_mantissa", param2.str(), inmap2.str(), outmap2.str());

		addComment("Adjust exponent");
		vhdl << declare(target->adderDelay(RegSize+wES), "sf_add", RegSize+wES) << " <= larger_sf + ovf_frac - lzCount;" << endl;

	//=========================================================================|
		addFullComment("Compute Regime and Exponent value");
	// ========================================================================|

		vhdl << declare(target->adderDelay(FracSize), "nzero", 1, false) << " <= '0' when frac_GR = \"" << string(FracSize+2, '0') << "\" else '1';"<<endl;

		addComment("Unpack scaling factors");
		if (wES > 0){
			vhdl << declare(.0, "FinalExp", wES) << " <= sf_add" << range(wES - 1, 0) << ";" << endl;
		}
		vhdl << declare(.0, "RegimeAns_tmp", RegSize) << " <= sf_add" << range(RegSize + wES - 1, wES) << ";" << endl;
		
		vhdl << declare(.0, "reg_sign", 1, false) << " <= RegimeAns_tmp(RegimeAns_tmp'high);" << endl;
		addComment("Get Regime's absolute value");

		if(wES>0){
			vhdl << "with reg_sign select " << declare(target->logicDelay() + target->adderDelay(RegSize), "FinalRegime", RegSize) << " <=" << endl <<
				tab << "(NOT RegimeAns_tmp) + 1 when '1'," << endl <<
				tab << "RegimeAns_tmp when '0'," << endl <<
				tab << "\"" << string(RegSize, '-') << "\" when others;" << endl;
		}
		else{
			addComment("Check for Regime overflow");

			std::string maxReg = toBinary(width-2, RegSize);

			vhdl << "with reg_sign select " << declare(target->logicDelay() + target->adderDelay(RegSize), "aux_FR", RegSize) << " <=" << endl <<
				tab << "(NOT RegimeAns_tmp)+1 when '1'," << endl <<
				tab << "RegimeAns_tmp when '0'," << endl <<
				tab << "\"" << string(RegSize, '-') << "\" when others;" << endl;

			/*
			vhdl << "with aux_FR" << range(RegSize-2,0) << " select " << declare(target->logicDelay() + target->adderDelay(RegSize), "FinalRegime", RegSize) << " <=" << endl <<
				tab << "'0' & \"" << string(RegSize-2, '1') << "\" & '0' when " << og(RegSize-1) << "," << endl <<
				tab << "aux_FR when others;" << endl;
			*/
			vhdl << declare(target->adderDelay(RegSize), "reg_ovf", 1, false) << " <= '1' when aux_FR > " << maxReg << " else '0';" << endl;
			vhdl << "with reg_ovf select " << declare(target->logicDelay(), "FinalRegime", RegSize) << " <=" << endl <<
				tab << maxReg << " when '1'," << endl <<
				tab << "aux_FR when '0'," << endl <<
				tab << "\"" << string(RegSize, '-') << "\" when others;" << endl;
		}
	
	//=========================================================================|
		addFullComment("Packing Stage 1");
	// ========================================================================|

	/*	vhdl << declare(.0, "tmp1", 2 + wES + FracSize + 2) << " <= nzero & '0' ";
		if (wES > 0){
			vhdl << "& FinalExp ";
		}
		vhdl << "& normFrac" << range(FracSize, 0) << "& S_tmp;" << endl;
		vhdl << declare(.0, "tmp2", 2 + wES + FracSize + 2) << " <= '0' & nzero ";
		if (wES > 0){
			vhdl << "& FinalExp ";
		}
		vhdl << "& normFrac" << range(FracSize, 0) << " & S_tmp;" << endl;

		vhdl << declare(target->adderDelay(RegSize-1), "shift_neg", RegSize-1) << " <= FinalRegime" << range(RegSize-2, 0) << " - 1;" << endl;
		vhdl << declare(target->adderDelay(RegSize-1), "shift_pos", RegSize-1) << " <= FinalRegime" << range(RegSize-2, 0) << " ;" << endl;

		vhdl << "with reg_sign select " << declare(target->logicDelay(), "input_shifter", 2 + wES + FracSize + 2) << "<=" << endl <<
			tab << "tmp2 when '1'," << endl <<
			tab << "tmp1 when '0'," << endl <<
			tab << "\"" << string(2 + wES + FracSize + 2, '-') << "\" when others;" << endl;
		vhdl << "with reg_sign select " << declare(target->logicDelay(), "shift_offset", RegSize-1) << " <=" << endl <<
			tab << "shift_neg when '1'," << endl <<
			tab << "shift_pos when '0'," << endl <<
			tab << "\"" << string(RegSize-1, '-') << "\" when others;" << endl;
		vhdl << declare(.0, "pad", 1, false) << "<= input_shifter(input_shifter'high);" << endl;
*/


		vhdl << "with reg_sign select " << declare(target->logicDelay(), "input_shifter", 2 + wES + FracSize + 2) << "<=" << endl <<
			tab << "'0' & nzero ";
			if (wES > 0){
				vhdl << tab << "& FinalExp ";
			}
				vhdl << tab << "& normFrac" << range(FracSize, 0) << " & S_tmp when '1'," << endl <<
			tab << "nzero & '0' ";
			if (wES > 0){
				vhdl << tab << "& FinalExp ";
			}
				vhdl << tab << "& normFrac" << range(FracSize, 0) << " & S_tmp when '0'," << endl <<
			tab << "\"" << string(2 + wES + FracSize + 2, '-') << "\" when others;" << endl;

		
		vhdl << "with reg_sign select " << declare(target->adderDelay(RegSize-1), "shift_offset", RegSize-1) << " <=" << endl <<
			tab << "FinalRegime" << range(RegSize-2, 0) << " - 1 when '1'," << endl <<
			tab << "FinalRegime" << range(RegSize-2, 0) << " when '0'," << endl <<
			tab << "\"" << string(RegSize-1, '-') << "\" when others;" << endl;
		
		vhdl << declare(.0, "pad", 1, false) << "<= input_shifter(input_shifter'high);" << endl;

		ostringstream param_shift, inmap_shift, outmap_shift;
		param_shift << "wIn=" << 2 + wES + FracSize + 2;
		param_shift << " maxShift=" << width-1;
		param_shift << " wOut=" << width + 1;	// width-1 + (G, R bits) PROBLEM: if(computeSticky) wOut =  wIn;
		param_shift << " dir=" << Shifter::Right;
		param_shift << " computeSticky=true";
		param_shift << " inputPadBit=true";
		
		inmap_shift << "X=>input_shifter,S=>shift_offset,padBit=>pad";
		
		outmap_shift << "R=>shifted_ans,Sticky=>S_bit_tmp";
		
		newInstance("Shifter", "right_signed_shifter", param_shift.str(), inmap_shift.str(), outmap_shift.str());
		int shift_size = getSignalByName("shifted_ans")->width(); // == 2 + wES + FracSize + 2

		vhdl << declare(.0, "tmp_ans", width-1) << " <= shifted_ans" << range(shift_size-1, shift_size-(width-1)) << ";" << endl;

	//=========================================================================|
		addFullComment("Packing Stage 2 - Unbiased Rounding");
	// ========================================================================|
		// Rounding implementation using L,G,R,S bits
		vhdl << declare(.0, "LSB", 1, false) << " <= shifted_ans" << of(shift_size-(width-1)) << ";" << endl;
		vhdl << declare(.0, "G_bit", 1, false) << " <= shifted_ans" << of(shift_size-(width-1)-1) << ";" << endl;
		vhdl << declare(.0, "R_bit", 1, false) << " <= shifted_ans" << of(shift_size-(width-1)-2) << ";" << endl;
		
//		vhdl << declare(target->adderDelay(shift_size-(width-1)-2), "S_bit", 1, false) << " <= S_bit_tmp when shifted_ans" << range(shift_size-(width-1)-3, 0) << " = \"" << string(shift_size-(width-1)-2, '0') << "\" else '1';" << endl;
		vhdl << declare(.0, "remain", shift_size-(width-1)-2) << " <= shifted_ans" << range(shift_size-(width-1)-3,0) << ";" << endl;
		//vhdl << declare(.0, "a2", shift_size-(width-1)-2) << " <= (others => '0');" << endl;
		vhdl << declare(target->adderDelay(shift_size-(width-1)-2), "S_bit", 1, false) << " <= S_bit_tmp when remain = \"" << string(shift_size-(width-1)-2, '0') << "\" else '1';" << endl;
		vhdl << declare(target->logicDelay(3), "round", 1, false) << " <= G_bit AND (LSB OR R_bit OR S_bit);" << endl;

		vhdl << "with larger_sign select " << declare(2* target->adderDelay(width-1) + target->logicDelay(), "result", width) << "<=" << endl <<
			tab << "'0' & (tmp_ans + round) when '0'," << endl <<
			tab << "'1' & ((NOT(tmp_ans + round))+1) when '1'," << endl <<
			tab << "\"" << string(width, '-') << "\" when others;" << endl;

		vhdl << "R <= '1' & \"" << string(width-1, '0') << "\" when inf = '1' else " << endl <<
			tab << zg(width) << " when nzero = '0' else" << endl <<
			tab << " result;" << endl;


		addFullComment("End of vhdl generation");

	};

	
	void PositAdd2::emulate(TestCase * tc) {
		// TODO: Change name of I/O signals to match with this code (Use X, Y, R)

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
		mpfr_add(r, x, y, GMP_RNDN);
		
		// Set outputs
		PositNumber posr(width, wES, r);
		mpz_class svR = posr.getSignalValue();
		tc->addExpectedOutput("R", svR);
		
		
		// clean up
		mpfr_clears(x, y, r, NULL);

	}


//	void PositAdd2::buildStandardTestCases(TestCaseList * tcl) {
		// please fill me with regression tests or corner case tests!
//	}


	OperatorPtr PositAdd2::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
		int width, wES;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width); 
    	UserInterface::parsePositiveInt(args, "wES", &wES);
    	return new PositAdd2(target, parentOp, width, wES);

	}

	void PositAdd2::registerFactory(){
		UserInterface::add("PositAdd2", // name
						"A posit adder with a single architecture.",
						"BasicPosit",
						"", //seeAlso
						"width(int): posit size in bits; \
                        				wES(int): exponent size in bits;",
						"",
						PositAdd2::parseArguments
						) ;
	}

}//namespace


