//
// Created by Viktor Schmidt.
//

#ifndef FLOPOCO_MONOTONEFUNCTION_H
#define FLOPOCO_MONOTONEFUNCTION_H

#include "ComparatorTable.hpp"
#include "FixMonotoneFunctionInterface.hpp"

namespace flopoco {
    class MonotoneFunctionComparator : public FixMonotoneFunctionInterface {

    public:
        MonotoneFunctionComparator(OperatorPtr parentOp, Target *target, string functionString_, int inputWidth, int outputWidth);

        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args);

        /** Factory register method */
        static void registerFactory();

        void build();
    };
}


#endif //FLOPOCO_MONOTONEFUNCTION_H
