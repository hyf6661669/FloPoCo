/*
  A Posit Decoder for FloPoCo

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
#include "PositDecoder.hpp"
#include <Operator.hpp>
#include <utils.hpp>

#include "ShiftersEtc/LZOC.hpp"
#include "ShiftersEtc/Shifters.hpp"

using namespace std;
namespace flopoco {
	
								
	PositDecoder::PositDecoder(Target* target, Operator* parentOp, int width_, int wES_):	Operator(parentOp, target), width(width_), wES(wES_) {
		srcFileName="PositDecoder";
		
		ostringstream name;

		name << "PositDecoder_" << width << "_" << wES;
		setNameWithFreqAndUID(name.str());
		setCopyrightString("Raul Murillo, Alberto A. del Barrio, Guillermo Botella, 2019");

		if (width < 3) {
			throw std::string("PositDecoder Constructor : width is too small, should be greater than two");
		}
		int freeWidth = width - 3;
		if (wES >= freeWidth) {
			//Avoid posits without even one bit of precision
			throw std::string("PositDecoder Constructor : invalid value of wES");
		}

		// SET UP THE IO SIGNALS

		// declaring inputs
		int sizeRegime = intlog2(width);
		int sizeFraction = width - wES-2;
		addInput  ( "Input"	, width);
		// declaring outputs
		addOutput ( "Sign"	);
		addOutput ( "Reg"	, sizeRegime);
		if(wES>0){
			addOutput ( "Exp"	, wES);
		}
		else{
			addOutput ( "Exp", 1);
		}
		addOutput ( "Frac"  , sizeFraction);
		addOutput ( "z" );
		addOutput ( "inf" );
		
		addFullComment("Start of vhdl generation");

		REPORT(INFO,"Declaration of PositDecoder \n");
		REPORT(DETAILED, "this operator has received two param1eters " << width << " and " << wES);
		REPORT(DEBUG,"debug of PositDecoder");

	//=========================================================================|
		addFullComment("Special Cases");
	// ========================================================================|
		
		vhdl << declare(target->adderDelay(width - 2), "nzero", 1, false) << " <= " <<
			"Input" << of(width - 2) << " when Input" << range(width - 3, 0) <<" = \"" << string(width -2, '0') << "\" else '1';" << endl;
		addComment("1 if Input is zero");
		vhdl << declare(target->logicDelay(2), "is_zero", 1, false) << " <= Input" << of(width - 1) << " NOR nzero;" << endl;	// 1 if Input is zero
		vhdl << "z <= is_zero;" << endl;
		addComment("1 if Input is infinity");
		vhdl << declare(target->logicDelay(2), "is_NAR", 1, false) << "<= Input" << of(width - 1) << " AND (NOT nzero);" << endl;	// 1 if Input is infinity
		vhdl << "inf <= is_NAR;" << endl;

	//=========================================================================|
		addFullComment("Extract Sign bit");
	// ========================================================================|
		vhdl << declare(.0, "my_sign", 1, false) << " <= Input" << of(width - 1) << ";" << endl;
		vhdl << "Sign <= my_sign;" << endl;

	//=========================================================================|
		addFullComment("2's Complement of Input");
	// ========================================================================|
		vhdl << declare(.0, "rep_sign", width - 1) << " <= (others => my_sign);" << endl;
		vhdl << declare(target->adderDelay(width - 1), "twos", width - 1) << " <= (rep_sign XOR Input" << range(width - 2, 0) << ") + my_sign;" << endl;
		vhdl << declare(.0, "rc", 1, false) << " <= twos" << of(width - 2) << ";" << endl;	// Regime check

	//=========================================================================|
		addFullComment("Count leading zeros of regime");
	// ========================================================================|
		// count the sequence of 0 bits terminating in a 1 bit - regime
		// zc ← Leading Zero Detector (inv)
		// we use the FloPoCo pipelined operator LZOC

		vhdl << declare(.0, "rep_rc", width - 1) << " <= (others => rc);" << endl;
		addComment("Invert 2's");
		vhdl << declare(target->logicDelay(1), "inv", width - 1) << " <= rep_rc XOR twos;" << endl;
	
		//vhdl << declare(.0, "zero_var", 1, false) << " <= '0';" << endl;
		ostringstream param1, inmap1, outmap1;
		param1 << "wIn=" << width - 1;
		param1 << " countType=" << 0;

		//inmap1 << "I=>inv,OZB=>zero_var";
		inmap1 << "I=>inv";
		
		outmap1 << "O=>zc";

		//newInstance("LZOC", "lzoc", param1.str(), inmap1.str(), outmap1.str());
		newInstance("LZOC", "lzc", param1.str(), inmap1.str(), outmap1.str());

	//=========================================================================|
		addFullComment("Shift out the regime");
	// ========================================================================|
		int zc_size = getSignalByName("zc")->width();
		vhdl << declare(target->adderDelay(zc_size), "zc_sub", zc_size) << " <= zc - 1;" << endl;
		ostringstream param2, inmap2, outmap2;
		param2 << "wIn=" << width - 1;
		param2 << " maxShift=" << width - 1;
		param2 << " dir=" << Shifter::Left;
    	//param2 << " inputPadBit=true";

		inmap2 << "X=>twos,S=>zc_sub";
		
		outmap2 << "R=>shifted_twos";

		newInstance("Shifter", "LeftShifterComponent", param2.str(), inmap2.str(), outmap2.str());

		vhdl << declare(.0, "tmp", width - 3) << " <= shifted_twos" << range(width - 4, 0) << ";"<<endl;

	//=========================================================================|
		addFullComment("Extract fraction and exponent");
	// ========================================================================|
		vhdl << "Frac <= nzero & tmp" << range(width - wES - 4, 0) << ";" << endl;
		if(wES>0){
			vhdl << "Exp <= tmp" << range(width - 4, width - wES - 3) << ";" << endl;
		}
		else{
			vhdl << "Exp <= \"0\";" << endl;
		}
		 

	//=========================================================================|
		addFullComment("Select regime");
	// ========================================================================|
		//vhdl << "Reg <= '0' & zc_sub when rc = '1' else " << endl <<
      	//	tab << "NOT('0' & zc) + 1;" << endl; //-zc

		vhdl << "with rc select " << declare(target->adderDelay(sizeRegime), "final_reg", sizeRegime) << "<= " << endl <<
      		tab << "'0' & zc_sub when '1'," << endl <<
			tab << "NOT('0' & zc) + 1 when '0'," << endl <<
			tab << "\"" << string(sizeRegime, '-') << "\" when others;" << endl;
		vhdl << "Reg <= final_reg;" << endl;

		addFullComment("End of vhdl generation");
	};

	
	void PositDecoder::emulate(TestCase * tc) {
	}


//	void PositDecoder::buildStandardTestCases(TestCaseList * tcl) {
		// please fill me with regression tests or corner case tests!
//	}


	OperatorPtr PositDecoder::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
		int width, wES;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width); 
    	UserInterface::parsePositiveInt(args, "wES", &wES);
    	return new PositDecoder(target, parentOp, width, wES);

	}
	
	void PositDecoder::registerFactory(){
		UserInterface::add("PositDecoder", // name
						"A posit decoder with a single architecture.",
						"BasicPosit",
						"", //seeAlso
						"width(int): posit size in bits; \
                        				wES(int): exponent size in bits;",
						"",
						PositDecoder::parseArguments
						) ;
	}

}//namespace
