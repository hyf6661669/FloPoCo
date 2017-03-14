#ifndef VariableColumnCompressor_HPP
#define VariableColumnCompressor_HPP

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
    class VariableColumnCompressor : public BasicCompressor
	{
	public:
		/** constructor **/
//        VariableColumnCompressor(Target * target, BasicCompressor *startCompressor, BasicCompressor *middleCompressor, BasicCompressor *endCompressor);
        VariableColumnCompressor(Target * target);

		/** destructor**/
		~VariableColumnCompressor();

        virtual void setWidth(int width) = 0;

    protected:
        int width;
	};

    std::ostream& operator<<(std::ostream& o, const VariableColumnCompressor& vcc);
}

#endif
