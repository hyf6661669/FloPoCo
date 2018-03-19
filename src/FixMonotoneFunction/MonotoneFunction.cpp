//
// Created by Viktor Schmidt.
//

#include "MonotoneFunction.hpp"

using namespace std;
namespace flopoco {
    MonotoneFunction::MonotoneFunction(OperatorPtr parentOp, Target* target, string functionString_, int inputWidth_, int outputWidth_)
            : FixMonotoneFunctionInterface(parentOp, target, functionString_, inputWidth_, outputWidth_) {
        srcFileName="FixMonotoneFunction";

        ostringstream name;
        name << "FixMonotoneFunction" << inputWidth << "_" << outputWidth;
        setName(name.str());

        REPORT(INFO,"Declaration of FixMonotoneFunction \n");

        // more detailed message
        REPORT(DETAILED, "this operator has received two parameters " << inputWidth << " and " << outputWidth);

        build();
    };


    mpz_class MonotoneFunction::calculateInverse(int y) {
        REPORT(DEBUG,"calculateInverse looking for x at f(x)=" << y);

        int ref = (int)pow(2, inputWidth - 1) -1;
        int nextOffset = (int)pow(2, inputWidth - 2);
        mpz_class lut_in(ref), lut_out;
        int sign = monotoneIncreasing ? 1 : -1;

        while(nextOffset != 0) {
            //lut_in = ;
            eval(lut_out, mpz_class(ref));

            int cmp = mpz_cmp_si(lut_out.get_mpz_t(), y);
            if(cmp > 0) {
                ref -= sign * nextOffset;
            }
            else {
                ref += sign * nextOffset;
            }

            //ref += mpz_cmp_si(lut_out.get_mpz_t(), y)

            nextOffset /= 2;

            REPORT(DEBUG,"nextoffset: " << nextOffset);
        }

        //lut_in = mpz_class(ref);

        eval(lut_out, mpz_class(ref));

        if(mpz_cmp_si(lut_out.get_mpz_t(), y)) {
            REPORT(DEBUG,"cmp non zero: " << mpz_cmp_si(lut_out.get_mpz_t(), y));
        }


//        ref += mpz_cmp_si(lut_out.get_mpz_t(), y);

        if(mpz_cmp_si(lut_out.get_mpz_t(), y) != 0) {
            ref += 1;
        }

        if(lut_in.get_si() > pow(2, inputWidth) - 1) {
            REPORT(DEBUG,"inverse too big: f(" << lut_in.get_str(2) << ")=" << y);
        }

        //lut_in = mpz_class(ref);

        return mpz_class(ref);
    }


    OperatorPtr MonotoneFunction::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
        string func;
        int inW, outW;
        UserInterface::parseString(args, "function", &func);
        UserInterface::parseInt(args, "inputWidth", &inW);
        UserInterface::parseInt(args, "outputWidth", &outW);
        return new MonotoneFunction(parentOp, target, func, inW, outW);
    }

    void MonotoneFunction::registerFactory(){
        UserInterface::add("MonotoneFunction", // name
                           "Generates a function.", // description, string
                           "Miscellaneous", // category, from the list defined in UserInterface.cpp
                           "",
                           "function(string)=x: Algorithm Type: normal, diff, lut;\
                           inputWidth(int)=16: Input bit count; \
                           outputWidth(int)=8: Output bit count",
                           "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developer manual in the doc/ directory of FloPoCo.",
                           MonotoneFunction::parseArguments
        ) ;
    }

    void MonotoneFunction::build() {
        string comparison = monotoneIncreasing ? ">=" : "<=";
        declare("output", outputWidth);
        vhdl << tab << "output(" << outputWidth - 1 << ") <= '1' when unsigned(i) " << comparison << " B\"" << calculateInverse(1 << (outputWidth - 1)).get_str(2) << "\" else '0';" << endl << endl;

        for(int i = 1; i < outputWidth; ++i) {
            std::vector<mpz_class> values = std::vector<mpz_class>();
            REPORT(DEBUG,"calculating LUT " << i);
            for(int j = 0; j < pow(2, i); ++j) {
                int v = (j << (outputWidth - i)) + (1 << (outputWidth - i - 1));
                values.push_back(calculateInverse(v));
            }

            ComparatorTable *ct = new ComparatorTable(this, getTarget(), i, inputWidth, values);
            addSubComponent(ct);

            string signal = declare(join("ref", i), inputWidth);

            ostringstream outputPort;
            outputPort << "output(" << outputWidth - 1 << " downto " << outputWidth - i << ")";
            this->inPortMap(ct, "X", outputPort.str());
            this->outPortMap(ct, "Y", signal);

            vhdl << this->instance(ct, join("ct", i));

            //vhdl << tab << "ct" << i << " : " << ct->getName() << " port map(X => output(" << outputWidth - 1 << " downto " << outputWidth - i << "), Y => " << signal << ", clk => clk, rst => rst);" << endl;
            vhdl << tab << "output(" << outputWidth - i - 1 << ") <= '1' when unsigned(i) " << comparison << " unsigned(" << signal << ") else '0';" << endl << endl;
        }

        vhdl << tab << "o <= output;" << endl;
    }


}//namespace