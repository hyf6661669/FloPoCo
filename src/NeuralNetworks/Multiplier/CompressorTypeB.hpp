#ifndef CompressorTypeB_HPP
#define CompressorTypeB_HPP

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
    class CompressorTypeB :public Operator
	{
	public:
		/** constructor **/
        CompressorTypeB(Target * target, int width_, int case3_);// case3: 3 => -A1; 4 => -A2; 5 => -A3; 6=> 0

		/** destructor**/
        ~CompressorTypeB();

        virtual void setWidth(int width_);
    private:
        bool useLastColumn;
        int width;
        int wOut;
        int case3;
       std::vector<int> outputs;
       std::vector<int> height;
	};
}
 
#endif
