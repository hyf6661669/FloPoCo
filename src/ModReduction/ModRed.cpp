#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "ModRed.hpp"

#include "PseudoCompressionStrategyOptILP.hpp"

using namespace std;

namespace flopoco {

    ModRed::ModRed(Operator *parentOp, Target *target_, int wIn_, int mod_) :
            Operator(parentOp, target_), wIn(wIn_), mod(mod_) {
        srcFileName = "ModRed";
        setCopyrightString("Andreas Boettcher");

        ostringstream name;
        name << "ModRed";
        setNameWithFreqAndUID(name.str());

        // the addition operators need the ieee_std_signed/unsigned libraries
        useNumericStd();

        //multiplierUid = parentOp->getNewUId();

        string xname = "X";

        // Set up the IO signals
        addInput(xname, wIn, true);
        int wOut = 1;
        for(; (1<<(wOut+1)) < mod; wOut++);
        addOutput("R", wOut, 2,true);

        BitHeap bitHeap(this, wIn);

        bitHeap.addSignal(xname, 0);
        //bitHeap.subtractSignal(xname, 0);

        CompressionStrategy *strategy = new PseudoCompressionptILP(&bitHeap, wIn, mod);
        bitHeap.startCompression(dynamic_cast<CompressionStrategy*>(strategy));

        vhdl << tab << "R" << " <= " << bitHeap.getSumName() << range(wOut, 0) << ";" << endl;

    }



    void ModRed::emulate (TestCase* tc)
    {
        mpz_class svX = tc->getInputValue("X");
        mpz_class svR;

        svR = svX%mod;
        tc->addExpectedOutput("R", svR);
    }

    void ModRed::buildStandardTestCases(TestCaseList* tcl)
    {
        TestCase *tc;

        mpz_class x, y;

        // 1
        x = mpz_class(1);
        tc = new TestCase(this);
        tc->addInput("X", x);
        emulate(tc);
        tcl->add(tc);

        // -1 mod
        x = (mpz_class(1) << wIn) -1;
        tc = new TestCase(this);
        tc->addInput("X", x);
        emulate(tc);
        tcl->add(tc);

        // The modulus of the two max negative value
        x = mpz_class(1) << (wIn -1);
        tc = new TestCase(this);
        tc->addInput("X", x);
        emulate(tc);
        tcl->add(tc);
    }

    OperatorPtr ModRed::parseArguments(OperatorPtr parentOp, Target *target, std::vector<std::string> &args) {
        int wIn_, mod_;

        UserInterface::parseStrictlyPositiveInt(args, "wIn", &wIn_);
        UserInterface::parseStrictlyPositiveInt(args, "mod", &mod_);


        return new ModRed(parentOp, target, wIn_, mod_);
    }

    void ModRed::registerFactory(){
        UserInterface::add("ModRed", // name
                           "A Implementation of a constant modulo reduction",
                           "BasicInteger", // category
                           "", // see also
                           "wIn(int): input size of the vector the constant modulo soul be calculated on;\
                        mod(int): constant modulus used of the calculation of teh remainder;", // This string will be parsed
                           "", // no particular extra doc needed
                           ModRed::parseArguments,
                           ModRed::unitTest
        ) ;
    }

    TestList ModRed::unitTest(int index)
    {
        // the static list of mandatory tests
        TestList testStateList;
        vector<pair<string,string>> paramList;

        list<int> wordSizes = {4,8,6,24,32,53,64,128,256};

        for(int modulus=1; modulus < 32; modulus+=2)
        {
            for (auto wordSizePair : wordSizes)
            {
                int wIn = wordSizePair;
                {
                    paramList.push_back(make_pair("wX", to_string(wIn)));
                    paramList.push_back(make_pair("mod", to_string(modulus)));
                    testStateList.push_back(paramList);
                    paramList.clear();
                }
            }
        }
        return testStateList;
    }



}