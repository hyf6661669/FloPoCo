/*
 * FixEMethodEvaluator.cpp
 *
 *  Created on: 7 Dec 2017
 *      Author: mistoan
 */

#include "FixEMethodEvaluator.hpp"

namespace flopoco {

	FixEMethodEvaluator::FixEMethodEvaluator(Target* target, size_t _radix, int _msbIn, int _lsbIn, int _msbOut, int _lsbOut,
		vector<mpfr_t> _coeffsP, vector<mpfr_t> _coeffsQ, map<string, double> inputDelays)
	: Operator(target), radix(_radix), n(_coeffsP.size()), m(coeffsQ.size()),
	  	  msbIn(_msbIn), lsbIn(_lsbIn), msbOut(_msbOut), lsbOut(_lsbOut),
		  maxDegree(n>m ? n : m)
	{
		ostringstream name;

		srcFileName = "FixEMethodEvaluator";
		name << "FixEMethodEvaluator_n" << n << "_m_" << m
				<< "_msbIn_" << vhdlize(msbIn) << "_lsbIn_" << vhdlize(lsbIn) << "_msbOut_" << vhdlize(msbOut) << "_lsbOut_" << vhdlize(lsbOut);
		setName(name.str());

		setCopyrightString("Matei Istoan, 2017");

		//safety checks and warnings
		if(n != m)
			REPORT(INFO, "WARNING: degree of numerator and of denominator are different! "
					+ "This will lead to a less efficient implementation.");

		wHatSize = -1;
		//W^, an estimation of W, made out of only a few MSBs of W
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

		//create a copy of the coefficients
		for(int i=0; i<_coeffsP.size(); i++)
		{
			mpfr_t tmpMpfr;

			mpfr_init2(tmpMpfr, _coeffsP[i]->_mpfr_prec);
			mpfr_set(tmpMpfr, _coeffsP[i], GMP_RNDN);
			coeffsP.push_back(tmpMpfr);
		}
		for(int i=_coeffsP.size(); i<maxDegree; i++)
		{
			mpfr_t tmpMpfr;

			mpfr_init2(tmpMpfr, 10000);
			mpfr_set_z(tmpMpfr, 0, GMP_RNDN);
			coeffsP.push_back(tmpMpfr);
		}
		for(int i=0; i<_coeffsQ.size(); i++)
		{
			mpfr_t tmpMpfr;

			mpfr_init2(tmpMpfr, _coeffsQ[i]->_mpfr_prec);
			mpfr_set(tmpMpfr, _coeffsQ[i], GMP_RNDN);
			coeffsQ.push_back(tmpMpfr);
		}
		for(int i=_coeffsQ.size(); i<maxDegree; i++)
		{
			mpfr_t tmpMpfr;

			mpfr_init2(tmpMpfr, 10000);
			mpfr_set_z(tmpMpfr, 0, GMP_RNDN);
			coeffsQ.push_back(tmpMpfr);
		}

		//compute the number of iterations needed
		nbIter = msbOut-lsbOut+1;
		if(radix > 2)
			nbIter = ceil(1.0*nbIter/intlog2(radix));
		//add an additional number of iterations to compensate for the errors
		g = intlog2(nbIter);
		nbIter += g;

		//add the inputs
		addFixInput("X", true, msbIn, lsbIn);
		//add the outputs
		addFixOutput("D", true, msbOut, lsbOut);

		//iteration 0
		addComment("iteration 0", ""+tab);
		for(size_t i=0; i<maxDegree; i++)
		{
			vhdl << tab << declareFixPoint(join("W_0_", i), true, msbOut, lsbOut-g) << " <= "
					<< signedFixPointNumber(coeffsP[i], msbIn, lsbIn, 0) << ";" << endl;
			vhdl << tab << declareFixPoint(join("D_0_", i), true, radix-1, 0) << " <= "
					<< zg(radix, 0) << ";" << endl;
		}

		//iterations 1 to nbIter-1
		for(size_t iter=1; iter<nbIter; iter++)
		{
			addComment(join("iteration ", iter), ""+tab);

			//first create element 0, which is simpler


			for(size_t i=0; i<maxDegree; i++)
			{

			}
		}
	}


	FixEMethodEvaluator::~FixEMethodEvaluator()
	{
		for(int i=0; i<coeffsP.size(); i++)
		{
			mpfr_clear(coeffsP[i]);
		}
		for(int i=0; i<coeffsQ.size(); i++)
		{
			mpfr_clear(coeffsQ[i]);
		}
	}

} /* namespace flopoco */
