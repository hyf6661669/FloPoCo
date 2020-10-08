//
// Created by annika on 31.08.20.
//

#include "ModuloMontgomery.hpp"

using namespace std;

namespace flopoco {
    ModuloMontgomery::ModuloMontgomery(OperatorPtr parentOp, Target* target, int wIn_, int modulo_, string method_)
    : Operator(parentOp,  target), wIn(wIn_), modulo(modulo_), method(method_) {
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

        // calculating R
        int modSize = getModuloMSB() + 1;
        R = (1 << modSize);

        // declaring inputs
        addInput ("X" , wIn, true);
        addConstant("M", "positive", modulo);
        addConstant("Rmod", "integer", R % modulo);
        // declaring output
        addOutput("S" , modSize);

        REPORT(INFO, "rmod: " << R % modulo);
        addFullComment("Start of vhdl generation");
        if(method == "modOp") {
            // calculating residue to use montgomery
            // conversion to signed
            vhdl << tab << declareFixPoint(
                    "XSg", false, modSize-1, 0) << tab << "<= UNSIGNED(STD_LOGIC_VECTOR(UNSIGNED(X) mod " << modulo << ")" << range(modSize-1, 0) << ");" << endl;
            vhdl << tab << declareFixPoint(
                    "RmodSg", false, modSize-1, 0) << tab << "<= UNSIGNED(TO_SIGNED(Rmod," << modSize << "));" << endl;
            // calculate and convert back
            vhdl << tab << declare(
                    "T_0", modSize*2*R, false) << tab << "<= (" << modSize*2*R-1 << " downto " << modSize*2 << " => '0') & STD_LOGIC_VECTOR(XSg * RmodSg);" << endl;
        } else if (method == "ex") {
            // classical modular multiplication
            // use shift because R is power of two
            vhdl << tab << declare(
                    "XRProd", wIn + modSize, false) << tab << "<= STD_LOGIC_VECTOR(SHIFT_LEFT(RESIZE(UNSIGNED(X)," << wIn + modSize << ")," << modSize << "));" << endl;
            // use division
            ostringstream divParams;
            divParams << "wIn=" << wIn + modSize << " d=" << modulo << " computeQuotient=false computeRemainder=true";
            newInstance("IntConstDiv", "modDiv", divParams.str(), "X=>XRProd", "R=>XRemain");
            vhdl << tab << declare(
                    "T_0", modSize*2*R, false) << tab << "<= (" << modSize*2*R-1 << " downto " << modSize << " => '0') & XRemain;" << endl;
        } else if (method == "redOnly") {
            if (wIn > modSize*2) {
                REPORT(INFO, "wIn too big for reduction only test");
            }

            vhdl << tab << declare(
                    "T_0", modSize*2*R, false) << tab << "<= (" << modSize*2*R-1 << " downto " << wIn << " => '0') & X;" << endl;
        } else {
            REPORT(INFO, "invalid method name " << method);
        }

        // montgomery reduction - automatically converts back
        for (int i = 0; i < modSize; ++i) {
            ostringstream tName;
            int index = i+1;
            tName << "T_" << index;
            vhdl << tab << declare(
                    tName.str(), modSize*2*R, false) << tab << "<= STD_LOGIC_VECTOR(UNSIGNED(T_" << i
                 << ") + SHIFT_LEFT(TO_UNSIGNED(M," << modSize+i << ")," << i << ")) when T_" << i << of(i) << " = '1' else T_" << i << ";" << endl;
        }
        vhdl << tab << declare(
                "TResTemp", modSize*2*R, false) << tab << "<= STD_LOGIC_VECTOR(SHIFT_RIGHT(UNSIGNED(T_" << modSize << ")," << modSize << "))"  << ";" << endl;
        vhdl << tab << declare(
                "TRes", modSize+1, false) << tab << " <= TResTemp" << range(modSize, 0) << ";" << endl;
        vhdl << tab << declare("TResMM", modSize+1, false) << "<= STD_LOGIC_VECTOR(UNSIGNED(TRes) - M)" << ";" << endl;
        vhdl << tab << "S <= TRes" << range(modSize-1, 0) << " when UNSIGNED(TRes) < M else TResMM" << range(modSize-1, 0) << ";" << endl;
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
            if ((number * i) % mod == 1) {
                return i;
            }
        }
    }

    void ModuloMontgomery::emulate(TestCase * tc) {
        mpz_class sx = tc->getInputValue("X");

        mpz_class sr;
        if (method == "redOnly") {
            int modSize = getModuloMSB() + 1;
            int rInverse = getModularInverse((1 << modSize), modulo);
            mpz_class sNormal = (sx * rInverse) % modulo;
            sr = sNormal % modulo;
        } else {
            sr = sx % modulo;
        }

        tc->addExpectedOutput("S",sr);
    }

    void ModuloMontgomery::buildStandardTestCases(TestCaseList * tcl) {
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

    TestList ModuloMontgomery::unitTest(int index)
    {
        // the static list of mandatory tests
        TestList testStateList;
        vector<pair<string,string>> paramList;
        string methods[] = {"ex","modOp","redOnly"};

        if(index==-1) 	{ // The unit tests
            for (int i = 0; i < 3; ++i) {
                for(int m=3; m<20; m+=2) {
                    int modSize = 0;
                    if (methods[i]=="redOnly") {
                        for (int w = 0; w < 31; ++w) {
                            if ((m & (1 << (w))) != 0) {
                                modSize = w;
                            }
                        }
                        modSize++;
                    }
                    for (int wIn=6; wIn<=12; wIn++) {
                        if (!(methods[i]=="redOnly" && wIn > modSize * 2)) {
                            paramList.push_back(make_pair("wIn",to_string(wIn)));
                            paramList.push_back(make_pair("modulo",to_string(m)));
                            paramList.push_back(make_pair("method", methods[i]));
                            paramList.push_back(make_pair("TestBench n=",to_string(100)));
                            testStateList.push_back(paramList);

                            paramList.clear();
                        }
                    }
                }
            }
        }
        return testStateList;
    }

    OperatorPtr ModuloMontgomery::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
        int wIn, modulo;
        string method;
        UserInterface::parseInt(args, "wIn", &wIn);
        UserInterface::parseInt(args, "modulo", &modulo);
        UserInterface::parseString(args, "method", &method);
        return new ModuloMontgomery(parentOp, target, wIn, modulo, method);
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
                            modulo(int): modulo; \
                            method(string)=ex: The method used for the montgomery algrithm."
                            " 3 available: 'ex' - for explicit conversion, 'modOp' - for Conversion using the vhdl mod operator"
                            " 'redOnly' - the input is directly used for the reduction without conversion in montgomery form",
                // More documentation for the HTML pages. If you want to link to your blog, it is here.
                           "See the developer manual in the doc/ directory of FloPoCo.",
                           ModuloMontgomery::parseArguments,
                           ModuloMontgomery::unitTest
        ) ;
    }
}
