//
// Created by Viktor Schmidt.
//

#include "MonotoneFunctionComparator.hpp"

using namespace std;
namespace flopoco {
    MonotoneFunctionComparator::MonotoneFunctionComparator(OperatorPtr parentOp, Target* target, string functionString_, int inputWidth_, int outputWidth_)
            : FixMonotoneFunctionInterface(parentOp, target, functionString_, inputWidth_, outputWidth_) {
        srcFileName="FixMonotoneFunctionComparator";

        ostringstream name;
        name << "FixMonotoneFunctionComparator" << inputWidth << "_" << outputWidth;
        setName(name.str());

        REPORT(INFO,"Declaration of FixMonotoneFunctionComparator \n");

        // more detailed message
        REPORT(DETAILED, "this operator has received two parameters " << inputWidth << " and " << outputWidth);

        build();
    };


    mpz_class MonotoneFunctionComparator::calculateInverse(int y) {
        REPORT(FULL,"calculateInverse looking for x at f(x)=" << y);

        long ref = (long)pow(2, inputWidth - 1);
        long nextOffset = (long)pow(2, inputWidth - 2);
        mpz_class lut_out;
        int sign = monotoneIncreasing ? 1 : -1;

        if(monotoneIncreasing) {
            --ref;
        }


        while(nextOffset != 0) {
            //lut_in = ;
            eval(lut_out, mpz_class(ref));

            //int cmp = mpz_cmp_si(lut_out.get_mpz_t(), y);

            if(lut_out.get_si() >= y) {
                ref -= sign * nextOffset;
            }
            else {
                ref += sign * nextOffset;
            }

            //ref += mpz_cmp_si(lut_out.get_mpz_t(), y)

            nextOffset /= 2;

            REPORT(FULL,"nextoffset: " << nextOffset);
        }

        //lut_in = mpz_class(ref);

        eval(lut_out, mpz_class(ref));


//        ref += mpz_cmp_si(lut_out.get_mpz_t(), y);

//        if(mpz_cmp_si(lut_out.get_mpz_t(), y) != 0) {
//            ref += 1;
//        }

        if(lut_out.get_si() < y) {
            ref += sign;
        }

//        if(lut_in.get_si() > pow(2, inputWidth) - 1) {
//            REPORT(DEBUG,"inverse too big: f(" << lut_in.get_str(2) << ")=" << y);
//        }

        //lut_in = mpz_class(ref);

        return mpz_class(ref);
    }


    OperatorPtr MonotoneFunctionComparator::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
        string func;
        int inW, outW;
        UserInterface::parseString(args, "function", &func);
        UserInterface::parseInt(args, "inputWidth", &inW);
        UserInterface::parseInt(args, "outputWidth", &outW);
        return new MonotoneFunctionComparator(parentOp, target, func, inW, outW);
    }

    void MonotoneFunctionComparator::registerFactory(){
        UserInterface::add("MonotoneFunctionComparator", // name
                           "Generates a function.", // description, string
                           "Miscellaneous", // category, from the list defined in UserInterface.cpp
                           "",
                           "function(string)=x: Algorithm Type: normal, diff, lut;\
                           inputWidth(int)=16: Input bit count; \
                           outputWidth(int)=8: Output bit count",
                           "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developer manual in the doc/ directory of FloPoCo.",
                           MonotoneFunctionComparator::parseArguments
        ) ;
    }

    void MonotoneFunctionComparator::build() {
        string comparison = monotoneIncreasing ? ">=" : "<=";
//        double delay = 0;

//        for(int x = 0; x < outputWidth; ++x) {
//            delay += getTarget()->adderDelay(inputWidth);
//            delay += getTarget()->tableDelay(x, inputWidth, true);
//        }


        //declare(getTarget()->localWireDelay(outputWidth), "output", outputWidth);
        string signal_comp_res = declare(getTarget()->adderDelay(inputWidth), "comp_res0", 1);
        //vhdl << tab << "output(" << outputWidth - 1 << ") <= '1' when unsigned(i) " << comparison << " B\"" << calculateInverse(1 << (outputWidth - 1)).get_str(2) << "\" else '0';" << endl << endl;
        vhdl << tab << signal_comp_res << "(0) <= '1' when unsigned(i) " << comparison << " B\"" << calculateInverse(1 << (outputWidth - 1)).get_str(2) << "\" else '0';" << endl << endl;

        for(int i = 1; i < outputWidth; ++i) {
            std::vector<mpz_class> values = std::vector<mpz_class>();
            REPORT(DEBUG,"calculating Table " << i);

            for(int j = 0; j < pow(2, i); ++j) {
                int v = (j << (outputWidth - i)) + (1 << (outputWidth - i - 1));
                values.push_back(calculateInverse(v));
            }

            ComparatorTable *ct = new ComparatorTable(this, getTarget(), i, inputWidth, values);
            addSubComponent(ct);

            string signal_ref = declare(getTarget()->tableDelay(i, inputWidth, true), join("ref", i), inputWidth);
            string signal_comp_res = declare(getTarget()->adderDelay(inputWidth) * i + getTarget()->tableDelay(i, inputWidth, true) * (i-1), join("comp_res", i), 1);
            string signal_table_in = declare(getTarget()->localWireDelay(i), join("table_in", i - 1), i);

            //ostringstream outputPort;
            //outputPort << "output(" << outputWidth - 1 << " downto " << outputWidth - i << ")";

            vhdl << tab << signal_table_in << " <= ";

            for(int k = 0; k < i; ++k) {
                vhdl << "comp_res" << k;

                if(k < i-1) {
                    vhdl << " & ";
                }
                else {
                    vhdl << ";" << endl << endl;
                }
            }

            this->inPortMap(ct, "X", signal_table_in);
            this->outPortMap(ct, "Y", signal_ref);

            vhdl << this->instance(ct, join("ct", i));

            //vhdl << tab << "ct" << i << " : " << ct->getName() << " port map(X => output(" << outputWidth - 1 << " downto " << outputWidth - i << "), Y => " << signal_ref << ", clk => clk, rst => rst);" << endl;
            //vhdl << tab << "output(" << outputWidth - i - 1 << ") <= '1' when unsigned(i) " << comparison << " unsigned(" << signal_ref << ") else '0';" << endl << endl;
            vhdl << tab << signal_comp_res << "(0) <= '1' when unsigned(i) " << comparison << " unsigned(" << signal_ref << ") else '0';" << endl << endl;
        }

        vhdl << tab << "o <= ";

        for(int k = 0; k < outputWidth; ++k) {
            vhdl << "comp_res" << k;

            if(k < outputWidth-1) {
                vhdl << " & ";
            }
            else {
                vhdl << ";" << endl << endl;
            }
        }

        //vhdl << tab << "o <= output;" << endl;
    }


}//namespace