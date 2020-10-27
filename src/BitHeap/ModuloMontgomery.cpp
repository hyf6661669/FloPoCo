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
        // declaring output
        addOutput("S" , modSize);

        addFullComment("Start of vhdl generation");

        if (method == "modOp" || method == "ex" || method == "redOnly") {
            // methods that convert the number before reduction to montgomery form
            if (method == "modOp") {
                addConstant("Rmod", "integer", R % modulo);
                REPORT(DEBUG, "rmod: " << R % modulo);
                // calculating residue to use montgomery
                // conversion to signed
                vhdl << tab << declareFixPoint(
                        "XSgTemp", false, wIn-1, 0) << tab << "<= UNSIGNED(STD_LOGIC_VECTOR(UNSIGNED(X) mod " << modulo << ")" << ");" << endl;
                vhdl << tab << declareFixPoint(
                        "XSg", false, modSize-1, 0) << tab << "<= XSgTemp" << range(modSize-1, 0) << ";" << endl;
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
                    THROWERROR("Parameter wIn is too big for reduction only test. Can only be twice the number of bits of the modulo.");
                }

                vhdl << tab << declare(
                        "T_0", modSize*2*R, false) << tab << "<= (" << modSize*2*R-1 << " downto " << wIn << " => '0') & X;" << endl;
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
        } else if (method == "multi") {
            // method that uses montgomery multiplication after the reduction to convert the result back
            if (wIn > modSize*2) {
                THROWERROR("Parameter wIn is too big for montgomery reduction. Can only be twice the number of bits of the modulo.");
            }

            int conversionMultiplicand = (1 << (2*modSize)) % modulo;
            addConstant("conMult", "natural", conversionMultiplicand);
            // montgomery reduction
            vhdl << tab << declare(
                    "T_0", modSize*2*R, false) << tab << "<= (" << modSize*2*R-1 << " downto " << wIn << " => '0') & X;" << endl;

            for (int i = 0; i < modSize; ++i) {
                ostringstream tName;
                int index = i+1;
                tName << "T_" << index;
                vhdl << tab << declare(getTarget()->adderDelay(modSize*2*R),
                        tName.str(), modSize*2*R, false) << tab << "<= STD_LOGIC_VECTOR(UNSIGNED(T_" << i
                     << ") + SHIFT_LEFT(TO_UNSIGNED(M," << modSize+i << ")," << i << ")) when T_" << i << of(i) << " = '1' else T_" << i << ";" << endl;
            }
            vhdl << tab << declare(
                    "TResTemp", modSize*2*R, false) << tab << "<= STD_LOGIC_VECTOR(SHIFT_RIGHT(UNSIGNED(T_" << modSize << ")," << modSize << "))"  << ";" << endl;
            vhdl << tab << declare(
                    "TRes", modSize+1, false) << tab << " <= TResTemp" << range(modSize, 0) << ";" << endl;
            vhdl << tab << declare(getTarget()->adderDelay(modSize+1),"TResMM", modSize+1, false) << "<= STD_LOGIC_VECTOR(UNSIGNED(TRes) - M)" << ";" << endl;
            vhdl << tab << declare(
                    "RedRes", modSize, false) << tab << " <= TRes" << range(modSize-1, 0) << " when UNSIGNED(TRes) < M else TResMM" << range(modSize-1, 0) << ";" << endl;

            // montgomery multiplication to convert back from montgomery form
            vhdl << tab << declare(
                    "R_0", modSize*2, false) << tab << "<= (" << modSize*2-1 << " downto 0 => '0');" << endl;
            for (int i = 0; i < modSize; ++i) {
                ostringstream rTmp1Name;
                int index = i+1;
                rTmp1Name << "Rtmp1_" << i;
                vhdl << tab << declare(getTarget()->adderDelay(modSize*2),
                        rTmp1Name.str(), modSize*2, false) << tab << "<= R_" << i << " when RedRes" << of(i) << " = '0' else STD_LOGIC_VECTOR(UNSIGNED(R_" << i << ") + conMult);" << endl;
                ostringstream rTmp2Name;
                rTmp2Name << "Rtmp2_" << i;
                vhdl << tab << declare(getTarget()->adderDelay(modSize*2),
                        rTmp2Name.str(), modSize*2, false) << tab << "<= Rtmp1_" << i << " when Rtmp1_" << i << of(0) << " = '0' else STD_LOGIC_VECTOR(UNSIGNED(Rtmp1_" << i << ") + M);" << endl;
                ostringstream rTmp3Name;
                rTmp3Name << "R_" << index;
                vhdl << tab << declare(
                        rTmp3Name.str(), modSize*2, false) << tab << "<= STD_LOGIC_VECTOR(SHIFT_RIGHT(UNSIGNED(Rtmp2_" << i << "),1));" << endl;
            }
            vhdl << tab << "S <= R_" << modSize << range(modSize-1,0) << " when UNSIGNED(R_" << modSize << ") < M";
            vhdl << tab << "else STD_LOGIC_VECTOR(UNSIGNED(R_" << modSize << range(modSize-1,0) << ") - M);" << endl;
            addFullComment("End of vhdl generation");

        } else if (method == "modred") {
            // montgomery reduction
            vhdl << tab << declare(
                    "T_0", wIn+1, false) << tab << "<= X & '0';" << endl;
            int length_red = 0;
            for (int i = 0; i < wIn; ++i) {
                ostringstream tName;
                tName << "T_" << i+1;
                vhdl << tab << declare(getTarget()->adderDelay(wIn+1-length_red),
                                       tName.str(), wIn+1-length_red, false) << tab << "<= STD_LOGIC_VECTOR(UNSIGNED( '0' & T_" << i
                     << "(" << wIn-length_red << " downto 1)) + TO_UNSIGNED(M, " << wIn+1-length_red << ")) when T_" << i << of(1) << " = '1' else '0' & T_" << i
                                                                                                             << "(" << wIn-length_red << " downto 1);" << endl;
                if(length_red+1 < wIn-modSize && 0 < i){
                    length_red++;
                }
            }
            vhdl << tab << declare(getTarget()->adderDelay(modSize+1),"TResMM", modSize, false) << "<= STD_LOGIC_VECTOR(UNSIGNED(T_" << wIn <<  "(" << ((wIn<modSize)?wIn:modSize) << " downto 1)) - M)" << ";" << endl;
            vhdl << tab << declare(
                    "RedRes", modSize, false) << tab << " <= T_" << wIn << range(modSize, 1) << " when UNSIGNED(T_" << wIn << range(wIn-length_red,1) <<  ") < M else TResMM" << ";" << endl;

            // montgomery multiplication to convert back from montgomery form

            //Constant to shift back from montgomery form
            mpz_class shiftConst = (((mpz_class)1 << 2*wIn) % modulo);
            vhdl << tab << declare(
                    "R_0", modSize+1, false) << tab << "<= (" << modSize << " downto 0 => '0');" << endl;
            for (int i = 0; i < wIn; ++i) {
                ostringstream rTmp1Name;
                rTmp1Name << "Rtmp1_" << i;
                int index = i + 1;
                mpz_class place_partial_mult = (shiftConst & ((mpz_class)1 << i)) >> i;
                if (place_partial_mult.get_si()){
                    vhdl << tab << declare(getTarget()->adderDelay(modSize+2),
                                           rTmp1Name.str(), modSize+2, false) << tab << "<= STD_LOGIC_VECTOR(UNSIGNED('0' & R_" << i << ") + UNSIGNED(\"00\" & RedRes));" << endl;
                }else {
                    vhdl << tab << declare(getTarget()->adderDelay(modSize+2),
                                           rTmp1Name.str(), modSize+2, false) << tab << "<= '0' & R_" << i << ";" << endl;
                }
                ostringstream rTmp2Name;
                rTmp2Name << "Rtmp2_" << i;
                vhdl << tab << declare(getTarget()->adderDelay(modSize+2),
                                       rTmp2Name.str(), modSize+2, false) << tab << "<= Rtmp1_" << i << " when Rtmp1_" << i << of(0) << " = '0' else STD_LOGIC_VECTOR(UNSIGNED(Rtmp1_" << i << ") + M);" << endl;
                ostringstream rTmp3Name;
                rTmp3Name << "R_" << index;
                vhdl << tab << declare(
                        rTmp3Name.str(), modSize+1, false) << tab << "<= STD_LOGIC_VECTOR(Rtmp2_" << i << range(modSize+1,1) << ");" << endl;
            }
            vhdl << tab << "S <= R_" << wIn << range(modSize-1,0) << " when UNSIGNED(R_" << wIn << ") < M";
            vhdl << tab << "else STD_LOGIC_VECTOR(UNSIGNED(R_" << wIn << range(modSize-1,0) << ") - M);" << endl;

            addFullComment("End of vhdl generation");

        } else if (method == "modred1") {
            // montgomery reduction
            vhdl << tab << declare(
                    "T_0", wIn+1, false) << tab << "<= X & '0';" << endl;

            int n = ((wIn-modSize)<modSize)?modSize:wIn-modSize;
            for (int i = 0; i < n; ++i) {
                ostringstream tName;
                tName << "T_" << i+1;
                vhdl << tab << declare(getTarget()->adderDelay(wIn+1),
                                       tName.str(), wIn+1, false) << tab << "<= STD_LOGIC_VECTOR(UNSIGNED( '0' & T_" << i
                     << "(" << wIn << " downto 1)) + TO_UNSIGNED(M, " << wIn+1 << ")) when T_" << i << of(1) << " = '1' else '0' & T_" << i
                     << "(" << wIn << " downto 1);" << endl;
            }
            vhdl << tab << declare(getTarget()->adderDelay(modSize+1),"TResMM", modSize, false) << "<= STD_LOGIC_VECTOR(UNSIGNED(T_" << n <<  "(" << ((wIn<modSize)?wIn:modSize) << " downto 1)) - M)" << ";" << endl;
            vhdl << tab << declare(
                    "RedRes", modSize, false) << tab << " <= T_" << n << range(modSize, 1) << " when UNSIGNED(T_" << n << range(wIn,1) <<  ") < M else TResMM" << ";" << endl;

            // montgomery multiplication to convert back from montgomery form

            //Constant to shift back from montgomery form
            mpz_class shiftConst = (((mpz_class)1 << 2*modSize) % modulo);
            vhdl << tab << declare(
                    "R_0", wIn, false) << tab << "<= (" << wIn-1 << " downto 0 => '0');" << endl;
            for (int i = 0; i < n; ++i) {
                ostringstream rTmp1Name;
                rTmp1Name << "Rtmp1_" << i;
                int index = i + 1;
                mpz_class place_partial_mult = (shiftConst & ((mpz_class)1 << i)) >> i;
                if (place_partial_mult.get_si()){
                    vhdl << tab << declare(getTarget()->adderDelay(wIn+1),
                                           rTmp1Name.str(), wIn+1, false) << tab << "<= STD_LOGIC_VECTOR(UNSIGNED('0' & R_" << i << ") + UNSIGNED('0' & RedRes));" << endl;
                }else {
                    vhdl << tab << declare(getTarget()->adderDelay(wIn+1),
                                           rTmp1Name.str(), wIn+1, false) << tab << "<= '0' & R_" << i << ";" << endl;
                }
                ostringstream rTmp2Name;
                rTmp2Name << "Rtmp2_" << i;
                vhdl << tab << declare(getTarget()->adderDelay(wIn+1),
                                       rTmp2Name.str(), wIn+1, false) << tab << "<= Rtmp1_" << i << " when Rtmp1_" << i << of(0) << " = '0' else STD_LOGIC_VECTOR(UNSIGNED(Rtmp1_" << i << ") + M);" << endl;
                ostringstream rTmp3Name;
                rTmp3Name << "R_" << index;
                vhdl << tab << declare(
                        rTmp3Name.str(), wIn, false) << tab << "<= STD_LOGIC_VECTOR(Rtmp2_" << i << range(wIn,1) << ");" << endl;
            }
            vhdl << tab << "S <= R_" << n  << range(modSize-1,0) << " when UNSIGNED(R_" << n  << ") < M";
            vhdl << tab << "else STD_LOGIC_VECTOR(UNSIGNED(R_" << n  << range(modSize-1,0) << ") - M);" << endl;

            addFullComment("End of vhdl generation");


        } else {
            REPORT(INFO, "invalid method name " << method);
        }
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
            // method multi
            for(int m=3; m<100; m+=2) {
                int modSize = 0;
                for (int w = 0; w < 31; ++w) {
                    if ((m & (1 << (w))) != 0) {
                        modSize = w;
                    }
                }
                modSize++;
                paramList.push_back(make_pair("wIn",to_string(modSize*2)));
                paramList.push_back(make_pair("modulo",to_string(m)));
                paramList.push_back(make_pair("method", "multi"));
                paramList.push_back(make_pair("TestBench n=",to_string(100)));
                testStateList.push_back(paramList);

                paramList.clear();
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
                            " 5 available: 'ex' - for explicit conversion, 'modOp' - for Conversion using the vhdl mod operator,"
                            " 'redOnly' - the input is directly used for the reduction without conversion in montgomery form,"
                            " 'multi' - using the montgomery reduction + montgomery multiplication to convert the result back"
                            " 'modred' - realizes modulo reduction for random sizes by montgomery reduction + montgomery multiplication to convert the result back",
                // More documentation for the HTML pages. If you want to link to your blog, it is here.
                           "See the developer manual in the doc/ directory of FloPoCo.",
                           ModuloMontgomery::parseArguments,
                           ModuloMontgomery::unitTest
        ) ;
    }
}
