//
// Created by annika on 31.08.20.
//

#include "../Operator.hpp"

#ifndef FLOPOCO_MODULOMONTGOMERY_H
#define FLOPOCO_MODULOMONTGOMERY_H

namespace flopoco {
    class ModuloMontgomery: public Operator {
        public:

            int wIn; // length of the input
            int modulo;
            string method;

        public:

            ModuloMontgomery(OperatorPtr parentOp, Target* target, int wIn = 3, int modulo = 1, string method = "ex");

            ~ModuloMontgomery() {};

            void emulate(TestCase * tc);

            void buildStandardTestCases(TestCaseList* tcl);

            static TestList unitTest(int index);

            /** Factory method that parses arguments and calls the constructor */
            static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);

            /** Factory register method */
            static void registerFactory();

        protected:

            int R;

            int getModuloMSB();

            int getModularInverse(int number, int mod);
    };
}//namespace

#endif //FLOPOCO_MODULOMONTGOMERY_H
