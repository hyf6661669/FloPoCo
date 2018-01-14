/*
 * GenericComputationUnit.cpp
 *
 *  Created on: 13 Jan 2018
 *      Author: mistoan
 */

#include "GenericComputationUnit.hpp"

namespace flopoco {

	GenericComputationUnit::GenericComputationUnit(Target* target, int _radix, int _maxDigit, int _index,
			Signal *_W, Signal *_X, Signal *_Di, string _qi, map<string, double> inputDelays)
	: Operator(target), radix(_radix), index(_index), maxDigit(_maxDigit),
	  msbW(_W->MSB()), lsbW(_W->LSB()),
	  msbX(_X->MSB()), lsbX(_X->LSB()),
	  msbD(_Di->MSB()), lsbD(_Di->LSB()),
	  qi(_qi)
	{
		ostringstream name;

		srcFileName = "GenericComputationUnit";
		name << "GenericComputationUnit_radix" << radix << "_index_" << index
				<< "_qi_" << std::setprecision(5) << qi << "_msbIn_" << vhdlize(msbW) << "_lsbIn_" << vhdlize(lsbW);
		setName(name.str());

		//safety checks
		if((radix != 2) && (radix != 4) && (radix != 8))
		{
			THROWERROR("GenericComputationUnit: radixes higher than 8 currently not supported!");
		}

		setCopyrightString("Matei Istoan, 2017");

		useNumericStd_Signed();

		//determine the MSB and the LSb for the internal computations
		msbInt = max(3, msbW, msbX, msbD);
		lsbInt = min(3, lsbW, lsbX, lsbD);

		//create the inputs and the output
		//	the inputs
		addFixInput("Wi", true, msbW, lsbW);
		addFixInput("D0", true, msbD, lsbD);
		addFixInput("Di", true, msbD, lsbD);
		addFixInput("Dip1", true, msbD, lsbD);
		addFixInput("X", true, msbX, lsbX);
		//	the inputs for D_{i+1}[j-1]*X
		for(int i=(-maxDigit); i<=maxDigit; i++)
			addFixInput(join("X_Mult_", vhdlize(i)), true, msbX+intlog2(maxDigit), lsbX);
		// the outputs
		addFixOutput("Wi_next", true, msbInt, lsbInt);

		//create the bitheap
		bitheap = new BitHeap(this, msbInt-lsbInt, false, join("Bitheap_"+name.str()+"_", getNewUId()));

		//add W_i[j-1] to the bitheap
		bitheap->addSignedBitVector(
									msbW-msbInt,					//weight
									"Wi",							//input signal name
									msbW-lsbW+1						//size
									);

		//create the multiplication D_0[j-1] * (-1)*q_i
		FixRealKCM *constMult = new FixRealKCM(
												this,				//parent operator
												"D0",				//input signal name
												true,				//signedness
												msbD,				//msbIn
												lsbD,				//lsbIn
												lsbInt,				//lsbOut
												"(-1)*"+qi,			//constant
												false,				//add round bit
												1.0					//target ulp error
												);
		//add the result of the multiplication to the bitheap
		constMult->addToBitHeap(bitheap, 0);

		//subtract D_i[j-1]
		bitheap->subtractSignedBitVector(
										msbD-msbInt,				//weight
										"Di",						//input signal name
										msbD-lsbD+1					//size
										);

		//create the multiplication D_{i+1}[j-1] * X
		//	not actual multiplication is done here, as we're multiplying by a digit
		//	instead, all the possible products are generated upstream
		//	and here we only have to choose which one to add
		vhdl << tab << declareFixPoint("Dip1_Mult_X", true, msbX+intlog2(maxDigit), lsbX) << " <= " << endl;
		for(int i=(-maxDigit); i<=maxDigit; i++)
		{
			mpz_class digitValue = i;

			//handle negative digits
			if(digitValue < 0)
				digitValue = mpz_class(1<<radix) + digitValue;

			vhdl << tab << tab << join("X_Mult_", vhdlize(i)) << " when Dip1="
					<< "\"" << unsignedBinary(digitValue, radix) << "\" else" << endl;
		}
		vhdl << tab << tab << "(others => '-');" << endl;
		//	add the result of the selection to the bitheap
		bitheap->addSignedBitVector(
									msbD-msbInt,					//weight
									"Dip1_Mult_X",					//input signal name
									msbD-lsbD+1						//size
									);

		//compress the bitheap
		bitheap->generateCompressorVHDL();

		// Retrieve the bits we want from the bit heap
		vhdl << declareFixPoint("sum", true, msbInt, lsbInt) << " <= " <<
				bitheap->getSumName() << range(msbInt-lsbInt, 0) << ";" << endl;

		vhdl << tab << "Wi_next <= sum;" << endl;
	}


	GenericComputationUnit::~GenericComputationUnit() {

	}


	void GenericComputationUnit::emulate(TestCase * tc)
	{

	}

	OperatorPtr GenericComputationUnit::parseArguments(Target *target, std::vector<std::string> &args) {
		int radix, index, maxDigit;
		int msbW, lsbW, msbX, lsbX, msbD, lsbD;
		string qi;

		UserInterface::parseInt(args, "radix", &radix);
		UserInterface::parseInt(args, "index", &index);
		UserInterface::parseInt(args, "maxDigit", &maxDigit);
		UserInterface::parseInt(args, "msbW", &msbW);
		UserInterface::parseInt(args, "lsbW", &lsbW);
		UserInterface::parseInt(args, "msbX", &msbX);
		UserInterface::parseInt(args, "lsbX", &lsbX);
		UserInterface::parseInt(args, "msbD", &msbD);
		UserInterface::parseInt(args, "lsbD", &lsbD);
		UserInterface::parseString(args, "qi", &qi);

		Signal *W  = new Signal("W", Signal::wire, true, msbW, lsbW);
		Signal *X  = new Signal("X", Signal::wire, true, msbX, lsbX);
		Signal *Di = new Signal("D", Signal::wire, true, msbD, lsbD);

		return new GenericComputationUnit(target, radix, index, maxDigit, W, X, Di, qi);
	}

	void GenericComputationUnit::registerFactory(){
		UserInterface::add("GenericComputationUnit", // name
				"Generic computation unit for the E-method.", //description
				"FunctionApproximation", // category
				"",
				"radix(int): the radix of the digit set being used;\
				 index(int): the index of the unit;\
				 maxDigit(int): the maximum digit in the used digit set;\
				 msbW(int): MSB of the W input signal;\
				 lsbW(int): LSB of the W input signal;\
				 msbX(int): MSB of the X input signal;\
				 lsbX(int): LSB of the X input signal;\
				 msbD(int): MSB of the D input signals;\
				 lsbD(int): LSB of the D input signals;\
				 q_i(string): the q_i constant, given in arbitrary-precision decimal, or as a Sollya expression, e.g \"log(2)\""
				"",
				"",
				GenericComputationUnit::parseArguments,
				GenericComputationUnit::unitTest
		) ;

	}

	TestList GenericComputationUnit::unitTest(int index)
	{
		// the static list of mandatory tests
		TestList testStateList;
		vector<pair<string,string>> paramList;

		return testStateList;
	}

} /* namespace flopoco */
