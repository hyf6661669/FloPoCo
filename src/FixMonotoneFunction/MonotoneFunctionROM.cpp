//
// Created by Viktor Schmidt.
//

#include "MonotoneFunctionROM.hpp"
// general c++ library for manipulating streams
#include <iostream>
#include <sstream>
#include <bitset>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"
#include "ComparatorTable.hpp"

using namespace std;
namespace flopoco {
    MonotoneFunctionROM::MonotoneFunctionROM(OperatorPtr parentOp, Target* target, string functionString_, int inputWidth_, int outputWidth_) :
            FixMonotoneFunctionInterface(parentOp, target, functionString_, inputWidth_, outputWidth_) {

        srcFileName="FixMonotoneFunctionROM";

        // definition of the name of the operator
        ostringstream name;
        name << "FixMonotoneFunctionROM" << inputWidth << "_" << outputWidth;
        setName(name.str());

        REPORT(INFO,"Declaration of FixMonotoneFunctionROM \n");

        REPORT(DETAILED, "this operator has received two parameters " << inputWidth << " and " << outputWidth);

        build();
    };


    OperatorPtr MonotoneFunctionROM::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
        string func;
        int inW, outW;
        UserInterface::parseString(args, "function", &func);
        UserInterface::parseInt(args, "inputWidth", &inW);
        UserInterface::parseInt(args, "outputWidth", &outW);
        return new MonotoneFunctionROM(parentOp, target, func, inW, outW);
    }

    void MonotoneFunctionROM::registerFactory(){
        UserInterface::add("MonotoneFunctionROM", // name
                           "Generates a LUT.", // description, string
                           "Miscellaneous", // category, from the list defined in UserInterface.cpp
                           "",
                           "function(string)=x: Algorithm Type: normal, diff, lut;\
                        inputWidth(int)=16: Input bit count; \
                        outputWidth(int)=8: Output bit count",
                           "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developer manual in the doc/ directory of FloPoCo.",
                           MonotoneFunctionROM::parseArguments
        ) ;
    }

    mpz_class MonotoneFunctionROM::function(int x) {
        mpz_class lut_in(x), lut_out;

        eval(lut_out, lut_in);

        return lut_out;
    }

    void MonotoneFunctionROM::build() {
//        addInput("X", inputWidth);
//        addOutput("Y", outputWidth);
//        string in = declare("inp", inputWidth);
//        string out = declare("outp", outputWidth);

        vector<mpz_class> values = vector<mpz_class>();

//        vhdl << tab << "with i select o <=" << endl;

        for(int x = 0; x < pow(2, inputWidth); ++x) {
            values.emplace_back(function(x));

//            vhdl << tab << tab << "\"" << unsignedBinary(function(x), outputWidth)
//                 << "\" when \"" << unsignedBinary(mpz_class(x), inputWidth) << "\"," << endl;

        }

//        vhdl << tab << tab << "\"";
//
//        for(int x = 0; x < outputWidth; ++x) {
//            vhdl << "-";
//        }

//        vhdl <<"\" when others;" << endl;
//        REPORT(DEBUG,"Last LUT value f(" << values.size() << ") = " << values.back().get_str(10) << " max y = " << y);

        ComparatorTable *ct = new ComparatorTable(this, getTarget(), inputWidth, outputWidth, values);

        //ct->inPortMap(this, in, "X");
        //ct->outPortMap(this, out, "Y");
        this->inPortMap(ct, "X", "i");
        this->outPortMap(ct, "Y", "o");
        addSubComponent(ct);

        vhdl << this->instance(ct, join("ct", inputWidth));
        //vhdl << tab << in << " <= i;" << endl;
        //vhdl << tab << "o <= " << out << ";" << endl;
    }

}//namespace