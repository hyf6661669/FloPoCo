//
// Created by Viktor Schmidt on 2/9/18.
//

#ifndef FLOPOCO_MONOTONEFUNCTIONDIFF_HPP
#define FLOPOCO_MONOTONEFUNCTIONDIFF_HPP

//#include "FixMonotoneFunction.hpp"
#include "FixMonotoneFunctionInterface.hpp"
#include "ComparatorTable.hpp"

namespace flopoco {
    class MonotoneFunctionDiff : public FixMonotoneFunctionInterface {

    public:
        MonotoneFunctionDiff(OperatorPtr parentOp, Target *target, string functionString_, int inputWidth, int outputWidth);

        //void emulate(TestCase *tc) override;

        //void buildStandardTestCases(TestCaseList *tcl) override;

        //void eval(mpz_class& r, mpz_class x) const;

        mpz_class calculateInverse(int y);
        void build();

        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args);
        static void registerFactory();
    };
}

#endif //FLOPOCO_MONOTONEFUNCTIONDIFF_HPP
