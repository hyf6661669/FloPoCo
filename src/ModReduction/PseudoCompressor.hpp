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
        PseudoCompressor(Operator *parentOp, Target *target, vector<int> _heights, vector<int> _outHeights);

        /** destructor**/
        ~PseudoCompressor();

        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args );
        static void registerFactory();

    protected:
        int _weight;
        int _modulus;
    };

    class BasicPseudoCompressor : public BasicCompressor
    {
    public:
        BasicPseudoCompressor(Operator* parentOp_, Target * target, vector<int> heights, vector<int> outHeights, int range_change, int _ones_vector_start=INT32_MAX);

        virtual Compressor* getCompressor(unsigned int middleLength);
    };
}
