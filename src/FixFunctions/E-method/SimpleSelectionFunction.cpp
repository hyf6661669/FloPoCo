/*
 * SimpleSelectionFunction.cpp
 *
 *  Created on: 8 Dec 2017
 *      Author: mistoan
 */

#include "SimpleSelectionFunction.hpp"

namespace flopoco {

		SimpleSelectionFunction::SimpleSelectionFunction(Target* target, int _radix, int _msbIn, int _lsbIn, map<string, double> inputDelays)
		: Operator(target), radix(_radix), msbIn(_msbIn), lsbIn(_lsbIn)
		{
			ostringstream name;

			name << "SimpleSelectionFunction_radix" << radix << "_msbIn_" << vhdlize(msbIn) << "_lsbIn_" << vhdlize(lsbIn);
			setName(name.str());
			//setNameWithFreqAndUID(name.str());
			setCopyrightString("Matei Istoan, 2017");

			//add the inputs
			addInput("W", msbIn-lsbIn+1);
			//add the outputs
			addOutput("D", radix);

			//create W^, an estimation of W, made out of only a few bits of W
			if(radix == 2)
			{
				//only the 4 top msbs are needed
				vhdl << tab << declare("WHat", 4) << " <= W" << range(msbIn-1, msbIn-5) << ";" << endl;
			}else
			{
				THROWERROR("SimpleSelectionFunction: radixes higher than 2 currently not supported");
			}

			vhdl << tab << "with WHat select D <= \n";
			for (mpz_class i = 0; i < (1 << radix); i++)
			{
				vhdl << tab << tab << "\"" << unsignedBinary((i+1)>>1, radix) << "\" when \""
						<< unsignedBinary(i, 4) << "\", \n";
			}
			vhdl << tab << tab << "\"" << std::string(radix, '-') << "\" when others;\n" << endl;

		}


		SimpleSelectionFunction::~SimpleSelectionFunction()
		{

		}

} /* namespace flopoco */
