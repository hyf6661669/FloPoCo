#ifndef BasicCompressor_HPP
#define BasicCompressor_HPP

#include "../Operator.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include "gmp.h"
#include "mpfr.h"
#include <vector>
#include <gmpxx.h>
#include <stdio.h>
#include <stdlib.h>
#include "../utils.hpp"

namespace flopoco
{

	/** The BasicCompressor class generates basic patterns for bit compressions
	 */


	class BasicCompressor:public Operator
	{
	public:
        vector<int> height; /** height of inputs, index 0 addresses the MSB column (!), e.g., a (1,5;3) GPC will have height[0]=1 and height[1]=5 **/
		int wOut; /** size of the output vector **/
		int param; /** computes the range of the output vector **/
		
		/* Uni KS start */
        vector<int> outputs; /** height of output bits, index 0 addresses the MSB column (!), for a GPC, each entry is 1, needed for compressors with several outputs per weight, e.g., 4:2 compressor **/
		double areaCost;        /** area cost of the compressor **/

		/** constructor **/
		BasicCompressor(Target * target);
		/* Uni KS stop */
		BasicCompressor(Target * target, vector<int> h);

		/** destructor**/
		~BasicCompressor();

		unsigned getColumnSize(int column);
		unsigned getNumberOfColumns();

        unsigned getOutputSize() const;


		/** test case generator  **/
		void emulate(TestCase * tc);

		// User-interface stuff
		/** Factory method */
		static OperatorPtr parseArguments(Target *target ,const vector<string> &args);

		static void registerFactory();
	};

    std::ostream& operator<<(std::ostream& o, const BasicCompressor& bc );

}

#endif
