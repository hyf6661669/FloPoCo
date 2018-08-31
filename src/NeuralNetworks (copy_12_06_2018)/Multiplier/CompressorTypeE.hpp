#ifndef CompressorTypeE_HPP
#define CompressorTypeE_HPP

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
    class CompressorTypeE :public Operator
	{
	public:
		/** constructor **/
        CompressorTypeE(Target * target, int width_, int case3_);// case3: 3 => -A1; 4 => -A2; 5 => -A3; 6=> 0

		/** destructor**/
        ~CompressorTypeE();

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
