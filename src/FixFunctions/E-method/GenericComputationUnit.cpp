/*
 * GenericComputationUnit.cpp
 *
 *  Created on: 13 Jan 2018
 *      Author: mistoan
 */

#include "GenericComputationUnit.hpp"

namespace flopoco {

	GenericComputationUnit::GenericComputationUnit(Target* target, int _radix, int _index,
			Signal *_W, Signal *_X, Signal *_Di, mpfr_t _qi, map<string, double> inputDelays)
	: Operator(target), radix(_radix), index(_index),
	  msbW(_W->MSB()), lsbW(_W->LSB()),
	  msbX(_X->MSB()), lsbX(_X->LSB()),
	  msbD(_Di->MSB()), lsbD(_Di->LSB())
	{
		ostringstream name;

		srcFileName = "GenericComputationUnit";
		name << "GenericComputationUnit_radix" << radix << "_index_" << index
				<< "_qi_" << std::setprecision(5) << mpfr_get_d(_qi, GMP_RNDN) << "_msbIn_" << vhdlize(msbW) << "_lsbIn_" << vhdlize(lsbW);
		setName(name.str());

		//safety checks
		if((radix != 2) && (radix != 4) && (radix != 8))
		{
			THROWERROR("GenericComputationUnit: radixes higher than 8 currently not supported!");
		}

		//make a copy of q_i
		mpfr_init2(qi, _qi->_mpfr_prec);
		mpfr_set(qi, _qi, GMP_RNDN);

		setCopyrightString("Matei Istoan, 2017");

		useNumericStd_Signed();

	}


	GenericComputationUnit::~GenericComputationUnit() {

	}


	void GenericComputationUnit::emulate(TestCase * tc)
	{

	}

	OperatorPtr GenericComputationUnit::parseArguments(Target *target, std::vector<std::string> &args) {
		return nullptr;
	}

	void GenericComputationUnit::registerFactory(){
		UserInterface::add("GenericComputationUnit", // name
				"Generic computation unit for the E-method.", //description
				"FunctionApproximation", // category
				"",
				"radix(int): the radix of the digit set being used;\
				 index(int): the index of the unit;\
				 msbW(int): MSB of the W input signal;\
				 lsbW(int): LSB of the W input signal;\
				 msbX(int): MSB of the X input signal;\
				 lsbX(int): LSB of the X input signal;\
				 msbD(int): MSB of the D input signals;\
				 lsbD(int): LSB of the D input signals;\
				 qi(string): the q_i constant, given in arbitrary-precision decimal, or as a Sollya expression, e.g \"log(2)\""
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
