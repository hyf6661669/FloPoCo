//
// Created by Annika Oeste on 24.10.20.
//

#include <math.h>
#include "ModuloBarrett.hpp"

using namespace std;
namespace flopoco {

    ModuloBarrett::ModuloBarrett(OperatorPtr parentOp, Target* target, int wIn_, int modulo_, int m_, int k_)
    : Operator(parentOp, target), wIn(wIn_), modulo(modulo_), m(m_), k(k_) {
        srcFileName="ModuloBarrett";

        // definition of the name of the operator
        ostringstream name;
        name << "ModuloBarrett" << wIn << "_" << modulo << "_" << m << "_" << k;
        setName(name.str()); // See also setNameWithFrequencyAndUID()
        // Copyright
        setCopyrightString("Annika Oeste, 2020");

        useNumericStd();

        REPORT(INFO,"Declaration of ModuloBarrett \n");
        REPORT(DETAILED, "this operator has received four parameters: wIn "
            << wIn << " modulo " << modulo << " m " << m << " k " << k);

        // declaring input
        addInput ("X" , wIn, true);
        // declaring output
        int modSize = floor(log2(modulo)+1);
        addOutput("Y" , modSize);

        addFullComment("Start of vhdl generation");

        // multiply with m
        ostringstream multMParams;
        multMParams << "wIn=" << wIn << " n=" << m;
        newInstance("IntConstMult", "XmMult", multMParams.str(), "X=>X", "R=>Xm");

        // only take msb without k lsb
        int xmSize = wIn + floor(log2(m)+1);
        int qSize = xmSize - k;

        vhdl << tab << declare(
                "Q", qSize, false) << tab << "<= Xm" << range(xmSize-1,k) << ";" << endl;
        // multiply with modulo
        ostringstream multModuloParams;
        multModuloParams << "wIn=" << qSize << " n=" << modulo;
        newInstance("IntConstMult", "QModuloMult", multModuloParams.str(), "X=>Q", "R=>QModulo");
        // subtract new result from X and put in Y
        int multSize = qSize + floor(log2(modulo)+1);
        vhdl << tab << declare(
                "YTmp", multSize, false) << tab << "<= STD_LOGIC_VECTOR(UNSIGNED(X) - UNSIGNED(QModulo));" << endl;
        vhdl << tab << "Y <= YTmp" << range(modSize-1, 0) << ";" << endl;
        addFullComment("End of vhdl generation");
    }

    void ModuloBarrett::emulate(TestCase * tc) {
        mpz_class sx = tc->getInputValue("X");
        mpz_class sr;

        sr = sx % modulo;
        tc->addExpectedOutput("Y",sr);
    }

    void ModuloBarrett::buildStandardTestCases(TestCaseList * tcl) {
        TestCase *tc;

        if (wIn >= 10) {
            tc = new TestCase(this);
            tc->addInput("X", mpz_class(634));
            emulate(tc);
            tcl->add(tc);
        }
        if (wIn >= 8) {
            tc = new TestCase(this);
            tc->addInput("X", mpz_class(222));
            emulate(tc);
            tcl->add(tc);
        }
        if (wIn >= 7) {
            tc = new TestCase(this);
            tc->addInput("X", mpz_class(120));
            emulate(tc);
            tcl->add(tc);
        }
        if (wIn >= 6) {
            tc = new TestCase(this);
            tc->addInput("X", mpz_class(33));
            emulate(tc);
            tcl->add(tc);

            tc = new TestCase(this);
            tc->addInput("X", mpz_class(40));
            emulate(tc);
            tcl->add(tc);
        }
        if (wIn >= 4) {
            tc = new TestCase(this);
            tc->addInput("X", mpz_class(15));
            emulate(tc);
            tcl->add(tc);
        }
    }

    TestList ModuloBarrett::unitTest(int index) {
        // the static list of mandatory tests
        TestList testStateList;
        vector<pair<string,string>> paramList;

        if(index==-1) 	{ // The unit tests
            for (int wIn=6; wIn<=12; wIn++) {
                for (int m = 2; m < 20; m++) {
                    paramList.push_back(make_pair("wIn",to_string(wIn)));
                    paramList.push_back(make_pair("modulo",to_string(m)));
                    paramList.push_back(make_pair("TestBench n=",to_string(100)));
                    testStateList.push_back(paramList);

                    paramList.clear();
                }
            }
        }

        return testStateList;
    }

    OperatorPtr ModuloBarrett::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
        int wIn, modulo, m, k;
        UserInterface::parseInt(args, "wIn", &wIn);
        UserInterface::parseInt(args, "modulo", &modulo);
        UserInterface::parseInt(args, "m", &m);
        UserInterface::parseInt(args, "k", &k);
        return new ModuloBarrett(parentOp, target, wIn, modulo, m, k);
    }

    void ModuloBarrett::registerFactory(){
        UserInterface::add("ModuloBarrett", // name
                           "An Operator to test the modulo calculation from P. Barrett", // description, string
                           "Miscellaneous", // category, from the list defined in UserInterface.cpp
                           "", //seeAlso
                // Now comes the parameter description string.
                // Respect its syntax because it will be used to generate the parser and the docs
                // Syntax is: a semicolon-separated list of parameterDescription;
                // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
                           "wIn(int)=16: A first parameter - the input size; \
                            modulo(int): modulo; \
                            m(int): constant m for the computation of S; \
                            k(int): constant k for the computation of S",
                // More documentation for the HTML pages. If you want to link to your blog, it is here.
                           "See the developer manual in the doc/ directory of FloPoCo.",
                           ModuloBarrett::parseArguments,
                           ModuloBarrett::unitTest
        ) ;
    }

}//namespace