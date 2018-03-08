//
// Created by Viktor Schmidt on 3/5/18.
//

#include "FixMonotoneFunction.hpp"

using namespace std;
namespace flopoco {
    OperatorPtr flopoco::FixMonotoneFunction::parseArguments(OperatorPtr parentOp, Target *target, vector <string> &args) {
        string param0;
        int param1, param2;
        UserInterface::parseString(args, "type", &param0);
        UserInterface::parseInt(args, "inputWidth", &param1);
        UserInterface::parseInt(args, "outputWidth", &param2);

        if(param0.compare("diff") == 0) {
            return new MonotoneFunctionDiff(target, param1, param2);
        }

        if(param0.compare("lut") == 0) {
            return new MonotoneFunctionLUT(target, param1, param2);
        }

        return new MonotoneFunction(target, param1, param2);
    }

    void FixMonotoneFunction::registerFactory() {
        UserInterface::add("FixMonotoneFunction", // name
                           "Generates a function.", // description, string
                           "Miscellaneous", // category, from the list defined in UserInterface.cpp
                           "", //seeAlso
                           "type(string)=normal: Algorithm Type: normal, diff, lut;\
                        inputWidth(int)=16: Input bit count; \
                        outputWidth(int)=8: Output bit count",
                           "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developer manual in the doc/ directory of FloPoCo.",
                           FixMonotoneFunction::parseArguments
        );
    }

    FixMonotoneFunction::FixMonotoneFunction(Operator *parentOp, Target *target) : Operator(parentOp, target) {}
}