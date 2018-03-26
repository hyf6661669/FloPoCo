//
// Created by Viktor Schmidt.
//

#ifndef FLOPOCO_MONOTONEFUNCTION_H
#define FLOPOCO_MONOTONEFUNCTION_H

#include "ComparatorTable.hpp"
#include "FixMonotoneFunctionInterface.hpp"

namespace flopoco {
    class MonotoneFunction : public FixMonotoneFunctionInterface {

    public:
        // definition of some function for the operator

        // constructor, defined there with two parameters (default value 0 for each)
        MonotoneFunction(OperatorPtr parentOp, Target *target, string functionString_, int inputWidth, int outputWidth);

        // destructor
        //~MonotoneFunction() {};

        /** Factory method that parses arguments and calls the constructor */
        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args);

        /** Factory register method */
        static void registerFactory();

        mpz_class calculateInverse(int y);

        void build();
    };
}


#endif //FLOPOCO_MONOTONEFUNCTION_H
