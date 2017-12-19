/*
 * SimpleSelectionFunction.cpp
 *
 *  Created on: 8 Dec 2017
 *      Author: Matei Istoan
 */

#include "SimpleSelectionFunction.hpp"

namespace flopoco {

		SimpleSelectionFunction::SimpleSelectionFunction(Target* target, int _radix, int maxDigit_, int _msbIn, int _lsbIn, map<string, double> inputDelays)
		: Operator(target), radix(_radix), maxDigit(maxDigit_), msbIn(_msbIn), lsbIn(_lsbIn)
		{
			ostringstream name;
			size_t wHatSize = -1;

			//safety checks
			if(maxDigit < 0)
				THROWERROR("maximum digit of the redundant digit set is negative!");
			if(maxDigit < (radix-1))
				REPORT(INFO, "WARNING: used digit set is not maximal!");
			if(intlog2(abs(maxDigit)) >= radix-1)
				THROWERROR("maximum digit larger than the maximum digit in the redundant digit set!");
			if(intlog2(abs(maxDigit)) >= (1<<msbIn))
				THROWERROR("maximum digit not representable on the given input format!");

			name << "SimpleSelectionFunction_radix" << radix << "_msbIn_" << vhdlize(msbIn) << "_lsbIn_" << vhdlize(lsbIn);
			setName(name.str());
			//setNameWithFreqAndUID(name.str());
			setCopyrightString("Matei Istoan, 2017");

			//add the inputs
			addInput("W", msbIn-lsbIn+1);
			//add the outputs
			addOutput("D", radix);

			//create W^, an estimation of W, made out of only a few MSBs of W
			if(radix == 2)
			{
				//only the 3 top msbs are needed
				vhdl << tab << declare("WHat", 3) << " <= W" << range(msbIn-1, msbIn-3) << ";" << endl;
				wHatSize = 3;
			}else if(radix == 4)
			{
				//only the 4 top msbs are needed
				vhdl << tab << declare("WHat", 4) << " <= W" << range(msbIn-1, msbIn-4) << ";" << endl;
				wHatSize = 4;
			}else if(radix == 8)
			{
				//only the 5 top msbs are needed
				vhdl << tab << declare("WHat", 5) << " <= W" << range(msbIn-1, msbIn-5) << ";" << endl;
				wHatSize = 5;
			}else
			{
				THROWERROR("SimpleSelectionFunction: radixes higher than 8 currently not supported!");
			}

			vhdl << tab << "with WHat select D <= \n";
			for(mpz_class i = mpz_class(-(1 << wHatSize)); i < (1 << wHatSize); i++)
			{
				int digitValue = (i+1) >> 1;
				if(digitValue < 0)
					digitValue = (1<<wHatSize) + digitValue;
				vhdl << tab << tab << "\"" << unsignedBinary(digitValue, radix) << "\" when \""
						<< unsignedBinary(i, wHatSize) << "\", \n";
			}
			vhdl << tab << tab << "\"" << std::string(radix, '-') << "\" when others;\n" << endl;

		}


		SimpleSelectionFunction::~SimpleSelectionFunction()
		{

		}

		void SimpleSelectionFunction::emulate(TestCase * tc)
		{
			// get the inputs from the TestCase
			mpz_class svW = tc->getInputValue("W");

			// manage signed digits
			mpz_class big1 = (mpz_class(1) << radix);
			mpz_class big1P = (mpz_class(1) << (radix-1));

			if(svW >= big1P)
				svW -= big1;

			// compute the multiple-precision output
			mpz_class svD = (svW + mpz_class(1)) >> 1;

			// manage two's complement at output
			if(svD < 0)
				svD += big1;

			// complete the TestCase with this expected output
			tc->addExpectedOutput("D", svD);
		}

		OperatorPtr SimpleSelectionFunction::parseArguments(Target *target, std::vector<std::string> &args) {
			int radix, maxDigit, msbIn, lsbIn;

			UserInterface::parseInt(args, "radix", &radix);
			UserInterface::parseInt(args, "maxDigit", &maxDigit);
			UserInterface::parseInt(args, "msbIn", &msbIn);
			UserInterface::parseInt(args, "lsbIn", &lsbIn);

			return new SimpleSelectionFunction(target, radix, maxDigit, msbIn, lsbIn);
		}

		void SimpleSelectionFunction::registerFactory(){
			UserInterface::add("BasicCompressor", // name
					"",
					"operator; floating point; floating-point adders", // categories
					"",
					"",
					SimpleSelectionFunction::parseArguments
			) ;

		}

} /* namespace flopoco */
