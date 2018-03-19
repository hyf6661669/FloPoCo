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

        if(type.compare("diff") == 0) {
            return new MonotoneFunctionDiff(parentOp, target, func, inW, outW);
        }

        if(type.compare("lut") == 0) {
            return new MonotoneFunctionLUT(parentOp, target, func, inW, outW);
        }

        return new MonotoneFunction(parentOp, target, func, inW, outW);
    }

    void FixMonotoneFunction::registerFactory() {
        UserInterface::add("FixMonotoneFunction", // name
                           "Generates a function.", // description, string
                           "Miscellaneous", // category, from the list defined in UserInterface.cpp
                           "", //seeAlso
                           "type(string)=normal: Algorithm Type: normal, diff, lut;\
                            function(string)=x: Algorithm Type: normal, diff, lut;\
                        inputWidth(int)=16: Input bit count; \
                        outputWidth(int)=8: Output bit count",
                           "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developer manual in the doc/ directory of FloPoCo.",
                           FixMonotoneFunction::parseArguments
        );
    }

    FixMonotoneFunction::FixMonotoneFunction(Operator *parentOp, Target *target) : Operator(parentOp, target) {}
}