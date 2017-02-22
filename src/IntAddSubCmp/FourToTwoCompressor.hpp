#ifndef FourToTwoCompressor_HPP
#define FourToTwoCompressor_HPP

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
#include "VariableColumnCompressor.hpp"


namespace flopoco
{
    class FourToTwoCompressor : public VariableColumnCompressor
	{
	public:
		/** constructor **/
        FourToTwoCompressor(Target * target, int width);

		/** destructor**/
        ~FourToTwoCompressor();

        virtual void setWidth(int width);
	};
}
 
#endif
