#ifndef TargetOptCompressor_HPP
#define TargetOptCompressor_HPP

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
#include "BasicCompressor.hpp"


namespace flopoco
{

	/** The TargetOptCompressor class provides a wrapper for target optimized compressors
	 */


	class TargetOptCompressor : public BasicCompressor
	{
	public:
		/** constructor **/
		TargetOptCompressor(Target * target, vector<int> h);

		/** destructor**/
		~TargetOptCompressor();

		vector<int> outputheight;

		unsigned getOutputSize();
	
	};
}
 
#endif
