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
    class MonotoneFunctionLUT : public Table {
    protected:
        /* operatorInfo is a user defined parameter (not a part of Operator class) for
           stocking information about the operator. The user is able to defined any number of parameter in this class, as soon as it does not affect Operator parameters undeliberatly*/
        int inputWidth;
        int outputWidth;

        sollya_obj_t fS;
    public:
        // definition of some function for the operator

        // constructor, defined there with two parameters (default value 0 for each)
        MonotoneFunctionLUT(OperatorPtr parentOp, Target *target, int inputWidth = 14, int outputWidth = 8);

        // destructor
        ~MonotoneFunctionLUT() {};

        mpz_class function(int x);


        // Below all the functions needed to test the operator
        /* the emulate function is used to simulate in software the operator
           in order to compare this result with those outputed by the vhdl opertator */
        void emulate(TestCase *tc);

        /* function used to create Standard testCase defined by the developper */
        void buildStandardTestCases(TestCaseList *tcl);


        /* function used to bias the (uniform by default) random test generator
           See FPExp.cpp for an example */
        // TestCase* buildRandomTestCase(int i);

        void eval(mpz_class& r, mpz_class x) const;
        void eval(mpfr_t r, mpfr_t x) const;

        /** Factory method that parses arguments and calls the constructor */
        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args);

        /** Factory register method */
        static void registerFactory();
    };
}

#endif //FLOPOCO_MONOTONELUT_H
