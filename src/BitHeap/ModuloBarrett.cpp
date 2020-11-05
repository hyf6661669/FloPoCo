//
// Created by Annika Oeste on 24.10.20.
//

#include <math.h>
#include "../utils.hpp"
#include "ModuloBarrett.hpp"

using namespace std;
namespace flopoco {

    ModuloBarrett::ModuloBarrett(OperatorPtr parentOp, Target* target, int wIn_, int modulo_, int m_, int k_)
    : Operator(parentOp, target), wIn(wIn_), modulo(modulo_), m(m_), k(k_) {
        srcFileName="ModuloBarrett";

        //calculation of coefficients m and k
        if(m < 0 || k < 0){
            double S_min = 1/(double)modulo;
            double S_max = ((1<<wIn)+modulo-1)/(double)((1<<wIn)+modulo-2)*1/(double)modulo;
            uint64_t m_min, m_max;
            k = wIn;
            while(true){
                m_min = (uint64_t)ceil(S_min*((uint64_t)1<<k));
                m_max = (uint64_t)floor(S_max*((uint64_t)1<<k));
                if(m_min <= m_max){
                   break;
                } else {
                    k++;
                }
            }
            m = (int)m_min;
            REPORT(DETAILED,"ModuloBarrett coefficients calculated: " << " m=" << m << " k=" << k);
        }

        // definition of the name of the operator
        ostringstream name;
        name << "ModuloBarrett" << wIn << "_" << modulo << "_" << m << "_" << k;
        setName(name.str()); // See also setNameWithFrequencyAndUID()
        // Copyright
        setCopyrightString("Annika Oeste, Andreas Boettcher, 2020");

        useNumericStd();

        REPORT(INFO,"Declaration of ModuloBarrett \n");
        REPORT(DETAILED, "this operator has received four parameters: wIn "
            << wIn << " modulo " << modulo << " m " << m << " k " << k);

        // declaring input
        addInput ("X" , wIn, true);
        // declaring output
        int modSize = floor(log2(modulo)+1);
        addOutput("Y" , modSize);

        // if wIn is smaller than the bit width of the modulo, no computation has to be done
        if (wIn < modSize) {
            addFullComment("Start of vhdl generation");
            vhdl << tab << "Y <= (" << modSize-1 << " downto " << wIn << " => '0') & X;" << endl;
            addFullComment("End of vhdl generation");
            return;
        }

        addFullComment("Start of vhdl generation");

        // multiply with m
        // extend X with 0 so that IntConstMultShiftAddOpt knows it is a positive number
        vhdl << tab << declare(
                "X0", wIn+1, false) << tab << "<= '0' & X;" << endl;
        ostringstream multMParams;
        multMParams << "wIn=" << wIn+1 << " constant=" << m;
        ostringstream outPortMapsXm;
        outPortMapsXm << "x_out0_c" << m << "=>Xm";
        newInstance("IntConstMultShiftAddOpt", "XmMult", multMParams.str(), "x_in0=>X0", outPortMapsXm.str());

        // only take msb without k lsb
        int xmSize = intlog2(m * ((mpz_class(1) << wIn)-1));
        int qSize = xmSize - k;

        // extend with 0 so that IntConstMultShiftAddOpt knows it is a positive number
        vhdl << tab << declare(
                "Q", qSize+1, false) << tab << "<= '0' & Xm" << range(xmSize-1,k) << ";" << endl;
        // multiply with modulo
        ostringstream multModuloParams;
        multModuloParams << "wIn=" << qSize+1 << " constant=" << modulo;
        ostringstream outPortMapsQModulo;
        outPortMapsQModulo << "x_out0_c" << modulo << "=>QModulo";
        newInstance("IntConstMultShiftAddOpt", "QModuloMult", multModuloParams.str(), "x_in0=>Q", outPortMapsQModulo.str());

        // subtract new result from X and put in Y
        int multSize = intlog2(modulo * ((mpz_class(1) << (qSize))-1));
        vhdl << tab << declare(
                "QModuloFitRange", multSize, false) << tab << "<= QModulo" << range(multSize-1,0) << ";" << endl;
        vhdl << tab << declare(
                "YTmp", multSize, false) << tab << "<= STD_LOGIC_VECTOR(UNSIGNED(X) - UNSIGNED(QModuloFitRange));" << endl;
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

        int wIns[] = {6,8,9,6,16,8,8};
        int modulos[] = {11,13,29,57,1009,20,16};
        int ms[] = {373,79,283,9,66511,205,1};
        int ks[] = {12,10,13,9,26,12,4};

        if(index==-1) 	{ // The unit tests
            for (int i = 0; i < 7; ++i) {
                paramList.push_back(make_pair("wIn",to_string(wIns[i])));
                paramList.push_back(make_pair("modulo",to_string(modulos[i])));
                paramList.push_back(make_pair("m",to_string(ms[i])));
                paramList.push_back(make_pair("k",to_string(ks[i])));
                paramList.push_back(make_pair("TestBench n=",to_string(100)));
                testStateList.push_back(paramList);

                paramList.clear();
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
                            m(int)=-1: constant m for the computation of S; \
                            k(int)=-1: constant k for the computation of S",
                // More documentation for the HTML pages. If you want to link to your blog, it is here.
                           "See the developer manual in the doc/ directory of FloPoCo.",
                           ModuloBarrett::parseArguments,
                           ModuloBarrett::unitTest
        ) ;
    }

}//namespace