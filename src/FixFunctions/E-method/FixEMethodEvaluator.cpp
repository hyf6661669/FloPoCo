/*
 * FixEMethodEvaluator.cpp
 *
 *  Created on: 7 Dec 2017
 *      Author: mistoan
 */

#include "FixEMethodEvaluator.hpp"

namespace flopoco {

	FixEMethodEvaluator::FixEMethodEvaluator(Target* target, size_t _radix, size_t _n, size_t _m, int _msbIn, int _lsbIn, int _msbOut, int _lsbOut,
		vector<mpfr_t> _coeffsP, vector<mpfr_t> _coeffsQ, map<string, double> inputDelays)
	: Operator(target), radix(_radix), n(_n), m(_m), msbIn(_msbIn), lsbIn(_lsbIn), msbOut(_msbOut), lsbOut(_lsbOut)
	{
		ostringstream name;

		srcFileName = "FixEMethodEvaluator";
		name << "FixEMethodEvaluator_n" << n << "_m_" << m
				<< "_msbIn_" << vhdlize(msbIn) << "_lsbIn_" << vhdlize(lsbIn) << "_msbOut_" << vhdlize(msbOut) << "_lsbOut_" << vhdlize(lsbOut);
		setName(name.str());

		setCopyrightString("Matei Istoan, 2017");

		//safety checks and warnings
		if(n != m)
			REPORT(INFO, "WARNING: degree of numerator and of denominator are different! This will lead to a less efficient implementation.");

		wHatSize = -1;
		//create W^, an estimation of W, made out of only a few MSBs of W
		if(radix == 2)
		{
			//only the 4 top MSBs are needed
			wHatSize = 4;
		}else if(radix == 4)
		{
			//only the 5 top MSBs are needed
			wHatSize = 5;
		}else if(radix == 8)
		{
			//only the 6 top MSBs are needed
			wHatSize = 6;
		}else
		{
			THROWERROR("FixEMethodEvaluator: radixes higher than 8 currently not supported!");
		}


	}


	FixEMethodEvaluator::~FixEMethodEvaluator()
	{

	}

} /* namespace flopoco */
