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

} /* namespace flopoco */
