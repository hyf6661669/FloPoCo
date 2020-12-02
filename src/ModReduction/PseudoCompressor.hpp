#pragma once

#include "Operator.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include "gmp.h"
#include "mpfr.h"
#include <vector>
#include <gmpxx.h>
#include <stdio.h>
#include <stdlib.h>
#include "utils.hpp"
#include "BitHeap/Compressor.hpp"


namespace flopoco
{
    class PseudoCompressor : public Compressor //: public VariableColumnCompressor
    {
    public:
        /** constructor **/
        PseudoCompressor(Operator* parentOp, Target* target, int weight, int modulus);

        /** destructor**/
        ~PseudoCompressor();

        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args );
        static void registerFactory();

    protected:
        int _weight;
        int _modulus;
    };
}
