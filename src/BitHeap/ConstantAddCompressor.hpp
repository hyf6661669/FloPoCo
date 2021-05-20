//
// Created by Annika Oeste on 06.05.21.
//
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

namespace flopoco {
    class ConstantAddCompressor: public Compressor {
    public:
        /** constructor **/
        ConstantAddCompressor(Operator* parentOp, Target* target, vector<int> _heights, int constant);

        /** destructor**/
        ~ConstantAddCompressor();

        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args );
        static void registerFactory();
    };

    class BasicConstantAddCompressor : public BasicCompressor
    {
    public:
        BasicConstantAddCompressor(Operator* parentOp_, Target * target, vector<int> _heights, int constant);

        virtual Compressor* getCompressor(unsigned int middleLength);

    protected:
        int _constant;
    };
}