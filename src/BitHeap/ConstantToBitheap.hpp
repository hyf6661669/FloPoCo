//
// Created by Annika Oeste on 27.07.21.
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
    class ConstantToBitheap: public Compressor {
    public:
        /** constructor **/
        ConstantToBitheap(Operator* parentOp, Target* target, int constant);

        /** destructor**/
        ~ConstantToBitheap();

        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args );
        static void registerFactory();
    };

    class BasicConstantToBitheap : public BasicCompressor
    {
    public:
        BasicConstantToBitheap(Operator* parentOp_, Target * target, int constant);

        virtual Compressor* getCompressor(unsigned int middleLength);

    protected:
        int _constant;
    };
}
