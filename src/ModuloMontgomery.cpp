//
// Created by annika on 31.08.20.
//

#include "ModuloMontgomery.hpp"

using namespace std;

namespace flopoco {
    ModuloMontgomery::ModuloMontgomery(OperatorPtr parentOp, Target* target, int wIn_, int modulo_) : Operator(parentOp,  target), wIn(wIn_), modulo(modulo_) {
        srcFileName="ModuloMontgomery";

        // definition of the name of the operator
        ostringstream name;
        name << "ModuloMontgomery" << wIn << "_" << modulo;
        setName(name.str()); // See also setNameWithFrequencyAndUID()
        // Copyright
        setCopyrightString("Annika Oeste, 2020");

        useNumericStd();

        REPORT(INFO,"Declaration of ModuloMontgomery \n");
        REPORT(DETAILED, "this operator has received two parameters: wIn " << wIn << " modulo " << modulo);

        // calculating R and m'
        R = (1 << wIn);
        mLine = -getModularInverse(R, modulo);

        // declaring inputs
        addInput ("X" , wIn, true);
        addConstant("M", "positive", modulo);
        addConstant("MLineInt", "integer", mLine);
        addConstant("Rmod", "integer", R % modulo);
        // declaring output
        int modSize = getModuloMSB() + 1;
        addOutput("S" , modSize*2);


        addFullComment("Start of vhdl generation");
        // calculating residue to use montgomery
        // conversion to signed
        vhdl << tab << declare(
                "MLine", modSize, false) << tab << "<= STD_LOGIC_VECTOR(TO_SIGNED(MLineInt," << modSize << "));" << endl;
        vhdl << tab << declareFixPoint(
                "XSg", true, modSize-1, 0) << tab << "<= SIGNED(STD_LOGIC_VECTOR(SIGNED(X) mod " << modulo << ")" << range(modSize-1, 0) << ");" << endl;
        vhdl << tab << declareFixPoint(
                "RmodSg", true, modSize-1, 0) << tab << "<= TO_SIGNED(Rmod," << modSize << ");" << endl;
        vhdl << tab << declareFixPoint(
                "test", true, modSize*2-1, 0) << tab << "<= XSg * RmodSg;" << endl;
        // calculate and convert back
        vhdl << tab << declare(
                "T", modSize*2, false) << tab << "<= STD_LOGIC_VECTOR(XSg * RmodSg);" << endl;
        // montgomery - automatically converts back

        for (int i = 0; i < modSize; ++i) {
            ostringstream uName;
            uName << "U_" << i;
            vhdl << tab << declare(
                    uName.str(), 1, false) << " <= MLine" << of(0) << " when T" << of(i) << " = '1' else '0';" << endl;
            vhdl << tab << "T <= STD_LOGIC_VECTOR(UNSIGNED(T) + SHIFT_LEFT(TO_UNSIGNED(M," << modSize << ")," << i << ")) when T" << of(i) << " = '1' else T;" << endl;
        }
        vhdl << tab << "T <= STD_LOGIC_VECTOR(SHIFT_RIGHT(UNSIGNED(T)," << modSize << "));" << endl;

        vhdl << tab << "S <= T when UNSIGNED(T) < M else STD_LOGIC_VECTOR(UNSIGNED(T) - M);" << endl;
        addFullComment("End of vhdl generation");
    }

    int ModuloMontgomery::getModuloMSB() {
        int mmsb = 0;
        for (int w = 0; w < 31; ++w) {
            if ((modulo & (1 << (w))) != 0) {
                mmsb = w;
            }
        }
        return mmsb;
    }

    int ModuloMontgomery::getModularInverse(int number, int mod) {
        for (int i = 0; i < mod; ++i) {
            if ((number * 1) % mod == 1) {
                return i;
            }
        }
    }

    void ModuloMontgomery::emulate(TestCase * tc) {
        mpz_class sx = tc->getInputValue("X");

        mpz_class sr;
        sr = sx % modulo;

        tc->addExpectedOutput("S",sr);
    }

    void ModuloMontgomery::buildStandardTestCases(TestCaseList * tcl) {
        TestCase *tc;

        if (wIn >= 8) {
            if (222 < modulo * R) {
                tc = new TestCase(this);
                tc->addInput("X", mpz_class(222));
                emulate(tc);
                tcl->add(tc);
            }
        } else if (wIn >= 7) {
            if (120 < modulo * R) {
                tc = new TestCase(this);
                tc->addInput("X", mpz_class(120));
                emulate(tc);
                tcl->add(tc);
            }
        } else if (wIn >= 6) {
            if (33 < modulo * R) {
                tc = new TestCase(this);
                tc->addInput("X", mpz_class(33));
                emulate(tc);
                tcl->add(tc);
            }

            if (40 < modulo * R) {
                tc = new TestCase(this);
                tc->addInput("X", mpz_class(40));
                emulate(tc);
                tcl->add(tc);
            }
        } else if (wIn >= 4) {
            if (15 < modulo * R) {
                tc = new TestCase(this);
                tc->addInput("X", mpz_class(15));
                emulate(tc);
                tcl->add(tc);
            }
        }
    }

    TestList ModuloMontgomery::unitTest(int index)
    {
        // the static list of mandatory tests
        TestList testStateList;
        vector<pair<string,string>> paramList;
        if(index==-1) 	{ // The unit tests
            for (int wIn=6; wIn<=12; wIn++) {
                for(int m=3; m<20; m+=2) {
                    paramList.push_back(make_pair("wIn",to_string(wIn)));
                    paramList.push_back(make_pair("modulo",to_string(m)));
                    paramList.push_back(make_pair("TestBench n=",to_string(1000)));
                    testStateList.push_back(paramList);

                    paramList.clear();
                }
            }
        }
        return testStateList;
    }

    OperatorPtr ModuloMontgomery::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
        int wIn, modulo;
        UserInterface::parseInt(args, "wIn", &wIn);
        UserInterface::parseInt(args, "modulo", &modulo);
        return new ModuloMontgomery(parentOp, target, wIn, modulo);
    }

    void ModuloMontgomery::registerFactory(){
        UserInterface::add("ModuloMontgomery", // name
                           "An Operator to test the modulo calculation from P. L. Montgomery", // description, string
                           "Miscellaneous", // category, from the list defined in UserInterface.cpp
                           "", //seeAlso
                // Now comes the parameter description string.
                // Respect its syntax because it will be used to generate the parser and the docs
                // Syntax is: a semicolon-separated list of parameterDescription;
                // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
                           "wIn(int)=16: A first parameter - the input size; \
                            modulo(int): modulo",
                // More documentation for the HTML pages. If you want to link to your blog, it is here.
                           "See the developer manual in the doc/ directory of FloPoCo.",
                           ModuloMontgomery::parseArguments,
                           ModuloMontgomery::unitTest
        ) ;
    }
}
