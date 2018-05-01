//
// Created by Viktor Schmidt on 3/5/18.
//

#include "FixMonotoneFunction.hpp"

using namespace std;
namespace flopoco {
    OperatorPtr flopoco::FixMonotoneFunction::parseArguments(OperatorPtr parentOp, Target *target, vector <string> &args) {
        string func, type;
        int inW, outW;
        UserInterface::parseString(args, "function", &func);
        UserInterface::parseString(args, "type", &type);
        UserInterface::parseInt(args, "inputWidth", &inW);
        UserInterface::parseInt(args, "outputWidth", &outW);

        transform(type.begin(), type.end(), type.begin(), ::tolower);

        if(type.compare("diff") == 0 || type.compare("difference") == 0) {
            return new MonotoneFunctionDiff(parentOp, target, func, inW, outW);
        }

        if(type.compare("rom") == 0) {
            return new MonotoneFunctionROM(parentOp, target, func, inW, outW);
        }

        return new MonotoneFunctionComparator(parentOp, target, func, inW, outW);
    }

    void FixMonotoneFunction::registerFactory() {
        UserInterface::add("FixMonotoneFunction", // name
                           "Generates a function.", // description, string
                           "Miscellaneous", // category, from the list defined in UserInterface.cpp
                           "", //seeAlso
                           "type(string)=normal: Algorithm Type: comp, diff, rom;\
                            function(string)=x: Algorithm Type: comp, diff, rom;\
                        inputWidth(int)=16: Input bit count; \
                        outputWidth(int)=8: Output bit count",
                           "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developer manual in the doc/ directory of FloPoCo.",
                           FixMonotoneFunction::parseArguments
        );
    }

    FixMonotoneFunction::FixMonotoneFunction(Operator *parentOp, Target *target) : Operator(parentOp, target) {}
}