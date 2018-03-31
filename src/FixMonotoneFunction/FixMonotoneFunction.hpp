//
// Created by Viktor Schmidt on 3/5/18.
//

#ifndef FLOPOCO_FIXMONOTONEFUNCTION_HPP
#define FLOPOCO_FIXMONOTONEFUNCTION_HPP

#include "Operator.hpp"
#include "utils.hpp"
#include <string>
#include <sstream>
#include "MonotoneFunctionComparator.hpp"
#include "MonotoneFunctionDiff.hpp"
#include "MonotoneFunctionROM.hpp"

namespace flopoco {
    class FixMonotoneFunction : public Operator {
    public:
        FixMonotoneFunction(Operator *parentOp, Target *target);

        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args);
        static void registerFactory();
    };
}


#endif //FLOPOCO_FIXMONOTONEFUNCTION_HPP
