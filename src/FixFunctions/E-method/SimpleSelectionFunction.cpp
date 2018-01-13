/*
 * SimpleSelectionFunction.cpp
 *
 *  Created on: 8 Dec 2017
 *      Author: Matei Istoan
 */

#include "SimpleSelectionFunction.hpp"

namespace flopoco {

		SimpleSelectionFunction::SimpleSelectionFunction(Target* target, int _radix, int maxDigit_, Signal *_W, map<string, double> inputDelays)
		: Operator(target), radix(_radix), maxDigit(maxDigit_), msbIn(_W->MSB()), lsbIn(_W->LSB())
		{
			ostringstream name;
			mpz_class iterLimit;

			srcFileName = "SimpleSelectionFunction";
			name << "SimpleSelectionFunction_radix" << radix << "_msbIn_" << vhdlize(msbIn) << "_lsbIn_" << vhdlize(lsbIn);
			setName(name.str());

			//safety checks
			if(maxDigit < 0)
				THROWERROR("maximum digit of the redundant digit set is negative!");
			if(maxDigit < (radix-1))
				REPORT(INFO, "WARNING: used digit set is not maximal!");
			if(intlog2(abs(maxDigit)) > radix-1)
				THROWERROR("maximum digit larger than the maximum digit in the redundant digit set!");
			if(intlog2(abs(maxDigit)) >= (1<<msbIn))
				THROWERROR("maximum digit not representable on the given input format!");

			setCopyrightString("Matei Istoan, 2017");

			useNumericStd_Signed();

			//get the format of the W^ signal
			getWHatFormat(radix, maxDigit, &msbWHat, &lsbWHat);
			wHatSize = msbWHat-lsbWHat+1;

			//safety check
			if((msbIn-lsbIn+1) < (int)wHatSize)
				THROWERROR("not enough digits on the given input format!");

			//add the inputs
			addFixInput("W", true, msbIn, lsbIn);
			//add the outputs
			addFixOutput("D", true, radix-1, 0);

			vhdl << tab << declareFixPoint("WHat", true, msbWHat, lsbWHat)
					<< " <= W" << range(msbIn-lsbIn, msbIn-lsbIn-wHatSize+1) << ";" << endl;

			vhdl << tab << "with WHat select D <= \n";
			if(((1<<wHatSize)-1) > (maxDigit+1))
				// maxDigit+1 so as to avoid possible overflows
				iterLimit = (maxDigit+1) << 1;
			else
				// maximum allowed number in this radix
				iterLimit = (1 << wHatSize) - 1;
			for(mpz_class i=-iterLimit; i<=iterLimit; i++)
			{
				mpz_class digitValue = (i+1) >> 1;
				mpz_class iAbs = i;

				//corner cases at the end of the intervals
				//	in order to keep the resulting digit in the allowed digit set
				if(digitValue > maxDigit)
					digitValue = maxDigit;
				else if(digitValue < -maxDigit)
					digitValue = -maxDigit;
				//handle negative digits at the output
				if(digitValue < 0)
					digitValue = mpz_class(1<<radix) + digitValue;
				//handle negative indices
				if(i < 0)
					iAbs = mpz_class(1<<wHatSize) + i;
				//create the corresponding output line
				vhdl << tab << tab << "\"" << unsignedBinary(digitValue, radix) << "\" when \""
						<< unsignedBinary(iAbs, wHatSize) << "\", \n";
			}
			vhdl << tab << tab << "\"" << std::string(radix, '-') << "\" when others;\n" << endl;

		}


		SimpleSelectionFunction::~SimpleSelectionFunction()
		{

		}


		void SimpleSelectionFunction::getWHatFormat(int _radix, int _maxDigit, int *_msb, int *_lsb)
		{
			size_t wHSize = 0;

			//W^ is an estimation of W, made out of only a few MSBs of W
			if(_radix == 2)
			{
				//only the 4 top MSBs are needed
				wHSize = 4;
			}else if(_radix == 4)
			{
				//only the 5 top MSBs are needed
				wHSize = 5;
			}else if(_radix == 8)
			{
				//only the 6 top MSBs are needed
				wHSize = 6;
			}else
			{
				cerr << "SimpleSelectionFunction: radixes higher than 8 currently not supported!";
				exit(1);
			}

			*(_lsb) = -1;
			*(_msb) = wHSize - 2;
		}


		void SimpleSelectionFunction::emulate(TestCase * tc)
		{
			// get the inputs from the TestCase
			mpz_class svW = tc->getInputValue("W");

			// manage signed digits
			mpz_class big1int = (mpz_class(1) << wHatSize);
			mpz_class big1intP = (mpz_class(1) << (wHatSize-1));
			mpz_class big1out = (mpz_class(1) << radix);

			if((msbIn-(int)wHatSize+1) > lsbIn)
			{
				if((msbIn-wHatSize) < 0)
					svW = svW >> abs((long int)((msbIn-wHatSize+1)-lsbIn));
				else
					svW = svW << abs((long int)((msbIn-wHatSize+1)-lsbIn));
			}

			if(svW >= big1intP)
				svW -= big1int;

			// compute the multiple-precision output
			mpz_class svD;
			// round down
			/*
			if(abs(svW) <= ((radix-1)<<1))
				svD = sgn(svW)*(abs(svW) + mpz_class(1)) >> 1;
			else
				svD = sgn(svW)*abs(svW) >> 1;
			*/
			// round to nearest
			if(abs(svW) <= ((radix-1)<<1))
				svD = (svW + mpz_class(1)) >> 1;
			else
				svD = sgn(svW)*abs(radix-1);

			// limit the output digit to the allowed digit set
			if(abs(svW) > maxDigit)
				svD = sgn(svW) * mpz_class(maxDigit);

			// manage two's complement at the output
			if(svD < 0)
				//svD += big1out;
				svD += big1int;
			//only use radix bits at the output
			svD = svD & ((1<<radix)-1);

			// complete the TestCase with this expected output
			tc->addExpectedOutput("D", svD);
		}

		OperatorPtr SimpleSelectionFunction::parseArguments(Target *target, std::vector<std::string> &args) {
			int radix, maxDigit, msbIn, lsbIn;

			UserInterface::parseInt(args, "radix", &radix);
			UserInterface::parseInt(args, "maxDigit", &maxDigit);
			UserInterface::parseInt(args, "msbIn", &msbIn);
			UserInterface::parseInt(args, "lsbIn", &lsbIn);

			Signal *W = new Signal("W", Signal::wire, true, msbIn, lsbIn);

			return new SimpleSelectionFunction(target, radix, maxDigit, W);
		}

		void SimpleSelectionFunction::registerFactory(){
			UserInterface::add("SimpleSelectionFunction", // name
					"Selection function for the E-method.", //description
					"FunctionApproximation", // category
					"",
					"radix(int): the radix of the digit set being used;\
						maxDigit(int): the maximum digit in the redundant digit set;\
						msbIn(int): MSB of the input;\
						lsbIn(int): LSB of the input"
					"",
					"",
					SimpleSelectionFunction::parseArguments,
					SimpleSelectionFunction::unitTest
			) ;

		}

		TestList SimpleSelectionFunction::unitTest(int index)
		{
			// the static list of mandatory tests
			TestList testStateList;
			vector<pair<string,string>> paramList;

			if(index==-1)
			{ // The unit tests
				paramList.push_back(make_pair("radix", to_string(2)));
				paramList.push_back(make_pair("maxDigit", to_string(1)));
				paramList.push_back(make_pair("msbIn", to_string(1)));
				paramList.push_back(make_pair("lsbIn", to_string(-1)));
				testStateList.push_back(paramList);
				paramList.clear();

				paramList.push_back(make_pair("radix", to_string(4)));
				paramList.push_back(make_pair("maxDigit", to_string(2)));
				paramList.push_back(make_pair("msbIn", to_string(2)));
				paramList.push_back(make_pair("lsbIn", to_string(-1)));
				testStateList.push_back(paramList);
				paramList.clear();

				paramList.push_back(make_pair("radix", to_string(4)));
				paramList.push_back(make_pair("maxDigit", to_string(3)));
				paramList.push_back(make_pair("msbIn", to_string(2)));
				paramList.push_back(make_pair("lsbIn", to_string(-1)));
				testStateList.push_back(paramList);
				paramList.clear();

				paramList.push_back(make_pair("radix", to_string(8)));
				paramList.push_back(make_pair("maxDigit", to_string(7)));
				paramList.push_back(make_pair("msbIn", to_string(3)));
				paramList.push_back(make_pair("lsbIn", to_string(-1)));
				testStateList.push_back(paramList);
				paramList.clear();
			}
			else
			{
				// finite number of random test computed out of index
			}

			return testStateList;
		}

} /* namespace flopoco */
