//
// Created by Viktor Schmidt.
//

#ifndef FLOPOCO_MONOTONEFUNCTION_H
#define FLOPOCO_MONOTONEFUNCTION_H

#include "ComparatorTable.hpp"
#include "FixMonotoneFunctionInterface.hpp"

namespace flopoco {
    class MonotoneFunction : public FixMonotoneFunctionInterface {
    protected:
        //int inputWidth;
        //int outputWidth;

        //bool monotoneIncreasing;

        //sollya_obj_t fS;
    public:
        // definition of some function for the operator

        // constructor, defined there with two parameters (default value 0 for each)
        MonotoneFunction(OperatorPtr parentOp, Target *target, string functionString_, int inputWidth, int outputWidth);

        // destructor
        ~MonotoneFunction() {};

        //void emulate(TestCase *tc) override;

        //void buildStandardTestCases(TestCaseList *tcl) override;

        //void eval(mpz_class& r, mpz_class x) const;

        /** Factory method that parses arguments and calls the constructor */
        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args);

        /** Factory register method */
        static void registerFactory();

        mpz_class calculateInverse(int y);

        void build();
    };
}


#endif //FLOPOCO_MONOTONEFUNCTION_H
