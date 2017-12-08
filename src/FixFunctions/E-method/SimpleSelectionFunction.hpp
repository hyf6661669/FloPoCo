/*
 * SimpleSelectionFunction.hpp
 *
 *  Created on: 8 Dec 2017
 *      Author: Matei Istoan
 */

#ifndef _SIMPLESELECTIONFUNCTION_HPP_
#define _SIMPLESELECTIONFUNCTION_HPP_

#include <vector>
#include <sstream>
#include <string>

#include <sollya.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "utils.hpp"

namespace flopoco {

	class SimpleSelectionFunction : public Operator
	{
	public:
		/**
		 * A simple constructor.
		 * Currently only implementing radix 2.
		 * Higher radixes (4, 8 etc.) to come.
		 * @param   radix          the radix being used
		 * @param   msbIn          MSB of the input
		 * @param   lsbIn          LSB of the input
		 */
		SimpleSelectionFunction(Target* target,
				int radix,
				int msbIn,
				int lsbIn,
				map<string, double> inputDelays = emptyDelayMap);

		/**
		 * Class destructor
		 */
		~SimpleSelectionFunction();

	private:
		int radix;                            /**< the radix being used */
		int msbIn;                            /**< MSB of the input */
		int lsbIn;                            /**< LSB of the input */

	};

} /* namespace flopoco */

#endif /* _SIMPLESELECTIONFUNCTION_HPP_ */
