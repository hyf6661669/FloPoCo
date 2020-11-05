//
// Created by Annika Oeste on 24.10.20.
//

#include "../Operator.hpp"

#ifndef FLOPOCO_MODULOBARRETT_H
#define FLOPOCO_MODULOBARRETT_H

namespace flopoco {
    class ModuloBarrett: public Operator {
    public:

        int wIn;
        int modulo;
        long long int m;
        int k;


    public:
        ModuloBarrett(OperatorPtr parentOp, Target* target, int wIn = 6, int modulo = 11, int m = 373, int k = 12);

        ~ModuloBarrett() {};

        void emulate(TestCase * tc);

        void buildStandardTestCases(TestCaseList* tcl);

        static TestList unitTest(int index);

        /** Factory method that parses arguments and calls the constructor */
        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);

        /** Factory register method */
        static void registerFactory();
    };
}//namespace

#endif //FLOPOCO_MODULOBARRETT_H
