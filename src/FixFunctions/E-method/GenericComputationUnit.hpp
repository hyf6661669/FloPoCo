/*
 * GenericComputationUnit.hpp
 *
 *  Created on: 13 Jan 2018
 *      Author: Matei Istoan
 */

#ifndef _GENERICCOMPUTATIONUNIT_HPP_
#define _GENERICCOMPUTATIONUNIT_HPP_

#include <vector>
#include <sstream>
#include <iostream>
#include <string>

#include <sollya.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "Signal.hpp"
#include "utils.hpp"

namespace flopoco {

	class GenericComputationUnit
	{
	public:
		/**
		 * A simple constructor.
		 *
		 * Performs the following computation:
		 *     w_i[j] = radix * (w_i[j-1] - d_0[j-1]*q_i - d_i[j-1] + d_{i+1}[j-1]*x)
		 * where
		 *     w_i[j] is the output of the computation unit
		 *     w_i[j-1], d_0[j-1], d_i[j-1], d_{i+1}[j-1], x are inputs
		 *     q_i is a constant, also an input
		 *     d_0, d_i and d_{i+1} are digits in the input radix
		 *
		 * Currently only implementing radix 2.
		 * Higher radixes (4, 8 etc.) to come.
		 * @param   radix          the radix being used
		 * @param   index          the index of the unit
		 * @param   W              the input signal W
		 * @param   X              the input signal X
		 * @param   Di             the input signal Di
		 * @param   qi             the coefficient q_i
		 */
		GenericComputationUnit(Target* target,
				int    radix,
				int    index,
				Signal *W,
				Signal *X,
				Signal *Di,
				mpfr_t qi,
				map<string, double> inputDelays = emptyDelayMap);

		/**
		 * Class destructor
		 */
		~GenericComputationUnit();


		/**
		 * Test case generator
		 */
		void emulate(TestCase * tc);

		// User-interface stuff
		/**
		 * Factory method
		 */
		static OperatorPtr parseArguments(Target *target, std::vector<std::string> &args);

		static void registerFactory();

		static TestList unitTest(int index);

	private:
		int radix;                            /**< the radix of the digit set being used */
		int index;                            /**< the index of the computation unit */

		int msbW;                             /**< MSB of the W signal */
		int lsbW;                             /**< LSB of the W signal */
		int msbX;                             /**< MSB of the X signal */
		int lsbX;                             /**< LSB of the X signal */
		int msbD;                             /**< MSB of the D signals */
		int lsbD;                             /**< LSB of the D signals */

		mpfr_t qi;                            /**< the q_i coefficient */
	};

} /* namespace flopoco */

#endif /* _GENERICCOMPUTATIONUNIT_HPP_ */
