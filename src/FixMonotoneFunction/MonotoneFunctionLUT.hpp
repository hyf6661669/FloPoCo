//
// Created by Viktor Schmidt.
//

#include "Table.hpp"
#include "utils.hpp"
#include <string>
#include <iomanip>
#include <sollya.h>
#include "FixMonotoneFunctionInterface.hpp"

#ifndef FLOPOCO_MONOTONELUT_H
#define FLOPOCO_MONOTONELUT_H

namespace flopoco {
    class MonotoneFunctionLUT : public FixMonotoneFunctionInterface {

    public:
        // definition of some function for the operator

        // constructor, defined there with two parameters (default value 0 for each)
        MonotoneFunctionLUT(OperatorPtr parentOp, Target *target, string functionString_, int inputWidth_, int outputWidth_);

        // destructor
        //~MonotoneFunctionLUT() {};

        mpz_class function(int x);


        /** Factory method that parses arguments and calls the constructor */
        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args);

        /** Factory register method */
        static void registerFactory();

        void build();
    };
}

#endif //FLOPOCO_MONOTONELUT_H
