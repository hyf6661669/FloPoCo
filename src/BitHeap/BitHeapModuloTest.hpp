//
// Created by Annika Oeste on 21.03.21.
//

#include "Operator.hpp"
#include "BitHeap/BitHeap.hpp"
#include "utils.hpp"

#ifndef FLOPOCO_BITHEAPMODULOTEST_H
#define FLOPOCO_BITHEAPMODULOTEST_H


namespace flopoco {

    // new operator class declaration
    class BitHeapModuloTest : public Operator {
    public:
        /* operatorInfo is a user defined parameter (not a part of Operator class) for
           stocking information about the operator. The user is able to defined any number of parameter in this class, as soon as it does not affect Operator parameters undeliberatly*/
        int wIn;
        int mod;
        int maxInput;

    public:
        // definition of some function for the operator

        // constructor, defined there with two parameters
        BitHeapModuloTest(Target* target, int wIn = 8, int mod = 11, int maxInput = -1, string mode = "default", string pseudoCompMode = "minRange");

        // destructor
        ~BitHeapModuloTest() {};

        int reqBitsForRange2Complement(int min, int max);

        // Below all the functions needed to test the operator
        /* the emulate function is used to simulate in software the operator
           in order to compare this result with those outputted by the vhdl operator */
        void emulate(TestCase * tc);

        /* function used to create Standard testCase defined by the developer */
        void buildStandardTestCases(TestCaseList* tcl);

        static TestList unitTest(int index);

        /* function used to bias the (uniform by default) random test generator
           See FPExp.cpp for an example */
        // TestCase* buildRandomTestCase(int i);

        /** Factory method that parses arguments and calls the constructor */
        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);

        /** Factory register method */
        static void registerFactory();
    };


}//namespace

#endif //FLOPOCO_BITHEAPMODULOTEST_H
