/*
  Posit decoder.

  Authors : Raul Murillo Montero, Alberto A. del Barrio, Guillermo Botella

  This file is part of the research article
  "Customized Posit Adders and Multipliers using the FloPoCo Core Generator"

  Initial software.
  Copyright © Complutense University of Madrid, Spain,
  2020.
  All rights reserved.
*/

#include <sstream>
#include <iostream>
#include "PositDecoder.hpp"
#include <Operator.hpp>
#include <utils.hpp>

//#include "../ShiftersEtc/LZOC.hpp"
//#include "../ShiftersEtc/Shifters.hpp"
//#include "ShiftersEtc/LZOC.hpp"
#include "ShiftersEtc/LZOCShifterSticky.hpp"

using namespace std;
namespace flopoco
{

	PositDecoder::PositDecoder(Target *target, Operator *parentOp, int width_, int wES_) : Operator(parentOp, target), width(width_), wES(wES_)
	{
		srcFileName = "PositDecoder";

		ostringstream name;

		name << "PositDecoder_" << width << "_" << wES;
		setNameWithFreqAndUID(name.str());
		setCopyrightString("Raul Murillo, 2020");

		if (width < 3)
		{
			throw std::string("PositDecoder Constructor : width is too small, should be greater than two");
		}
		int freeWidth = width - 3;
		if (wES >= freeWidth)
		{
			//Avoid posits without even one bit of precision
			throw std::string("PositDecoder Constructor : invalid value of wES");
		}

		// SET UP THE IO SIGNALS

		// declaring inputs
		int sizeRegime = intlog2(width - 1) + 1;
		int sizeFraction = width - wES - 2;
		addInput("Input", width);
		// declaring outputs
		addOutput("Sign");
		addOutput("Reg", sizeRegime);
		if (wES > 0)
		{
			addOutput("Exp", wES);
		}
		else
		{
			addOutput("Exp", 1);
		}
		addOutput("Frac", sizeFraction);
		addOutput("z");
		addOutput("inf");
		addOutput("Abs_in", width - 1);

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of PositDecoder \n");
		REPORT(DETAILED, "this operator has received two param1eters " << width << " and " << wES);
		REPORT(DEBUG, "debug of PositDecoder");

		//=========================================================================|
		addFullComment("Extract Sign bit");
		// ========================================================================|
		vhdl << declare(.0, "s", 1, false) << " <= Input" << of(width - 1) << ";" << endl;
		vhdl << "Sign <= s;" << endl;

		//=========================================================================|
		addFullComment("Special Cases");
		// ========================================================================|

		vhdl << declare(target->adderDelay(width - 2), "nzero", 1, false) << " <= "
			 << "Input" << of(width - 2) << " when Input" << range(width - 3, 0) << " = \"" << string(width - 2, '0') << "\" else '1';" << endl;
		addComment("1 if Input is zero");
		vhdl << declare(target->logicDelay(), "is_zero", 1, false) << " <= s NOR nzero;" << endl; // 1 if Input is zero
		vhdl << "z <= is_zero;" << endl;
		addComment("1 if Input is infinity");
		vhdl << declare(target->logicDelay(2), "is_NAR", 1, false) << "<= s AND (NOT nzero);" << endl; // 1 if Input is infinity
		vhdl << "inf <= is_NAR;" << endl;

		//=========================================================================|
		addFullComment("2's Complement of Input");
		// ========================================================================|
		vhdl << declare(.0, "rep_sign", width - 1) << " <= (others => s);" << endl;
		vhdl << declare(target->logicDelay() + target->adderDelay(width - 1), "twos", width - 1) << " <= (rep_sign XOR Input" << range(width - 2, 0) << ") + s;" << endl;
		vhdl << declare(.0, "rc", 1, false) << " <= twos(twos'high);" << endl; // Regime check

		//=========================================================================|
		addFullComment("Count leading zeros of regime & shift it out");
		// ========================================================================|
		// count the sequence of 0 bits terminating in a 1 bit - regime
		// zc ← Leading Zero Detector (inv)
		// we use the FloPoCo pipelined operator LZOC
		vhdl << declare(.0, "remainder", width - 2) << "<= twos" << range(width - 3, 0) << ";" << endl;

		ostringstream param, inmap, outmap;
		int wCount = intlog2(width - 2); //sizeRegime - 1;	// intlog2(width-1)

		param << "wIn=" << width - 2;
		param << " wOut=" << wES + sizeFraction;
		param << " wCount=" << wCount;

		inmap << "I=>remainder,OZb=>rc";

		outmap << "Count=>lzCount,O=>usefulBits";

		newInstance("LZOCShifterSticky", "lzoc", param.str(), inmap.str(), outmap.str());
		int wShifted = getSignalByName("usefulBits")->width();

		//=========================================================================|
		addFullComment("Extract fraction and exponent");
		// ========================================================================|
		vhdl << "Frac <= nzero & usefulBits" << range(wShifted - wES - 2, wShifted - wES - sizeFraction) << ";" << endl;
		if (wES > 0)
		{
			vhdl << "Exp <= usefulBits" << range(wShifted - 2, wShifted - wES - 1) << ";" << endl;
		}
		else
		{
			vhdl << "Exp <= \"0\";" << endl;
		}

		//=========================================================================|
		addFullComment("Select regime");
		// ========================================================================|
		//vhdl << "Reg <= '0' & zc_sub when rc = '1' else " << endl <<
		//	tab << "NOT('0' & zc) + 1;" << endl; //-zc

		// lzCount = #(sequence of bits with the same value)-1
		vhdl << "with rc select " << declare(target->logicDelay(), "final_reg", sizeRegime) << "<= " << endl
			 << tab << "\"" << string(sizeRegime - wCount, '0') << "\" & lzCount when '1'," << endl
			 <<
			//tab << "NOT('0' & (lzCount+1)) + 1 when '0'," << endl <<
			tab << "NOT(\"" << string(sizeRegime - wCount, '0') << "\" & lzCount)  when '0'," << endl
			 << // equivalent expression as commented line above
			tab << "\"" << string(sizeRegime, '-') << "\" when others;" << endl;
		vhdl << "Reg <= final_reg;" << endl;
		vhdl << "Abs_in <= twos;" << endl;

		addFullComment("End of vhdl generation");
	};

	void PositDecoder::emulate(TestCase *tc)
	{
	}

	//	void PositDecoder::buildStandardTestCases(TestCaseList * tcl) {
	// please fill me with regression tests or corner case tests!
	//	}

	OperatorPtr PositDecoder::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int width, wES;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width);
		UserInterface::parsePositiveInt(args, "wES", &wES);
		return new PositDecoder(target, parentOp, width, wES);
	}

	void PositDecoder::registerFactory()
	{
		UserInterface::add("PositDecoder", // name
						   "A posit decoder with a single architecture.",
						   "BasicPosit",
						   "", //seeAlso
						   "width(int): posit size in bits; \
                        		wES(int): exponent size in bits;",
						   "",
						   PositDecoder::parseArguments);
	}

} // namespace flopoco
