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
                           "function(string)=x: Sollya function;\
                           inputWidth(int)=16: Input bit count; \
                           outputWidth(int)=8: Output bit count",
                           "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developer manual in the doc/ directory of FloPoCo.",
                           MonotoneFunctionComparator::parseArguments
        ) ;
    }

    void MonotoneFunctionComparator::build() {
        string comparison = monotoneIncreasing ? ">=" : "<";
        mpz_class inputWidth_ref = mpz_class(1) << inputWidth;

        mpz_class inverse = calculateInverse(1 << (outputWidth - 1));
        mpz_class lut_out;

        eval(lut_out, inverse);

        if(!monotoneIncreasing && lut_out.get_si() == 1 << (outputWidth - 1)) {
            ++inverse;
        }

        vhdl << tab << declare(getTarget()->adderDelay(inputWidth), "comp_res0", 1)
             << "(0) <= '1' when unsigned(i) " << comparison
             << " B\"" << inverse.get_str(2)
             << "\" else '0';" << endl << endl;

        for(int i = 1; i < outputWidth; ++i) {
            std::vector<mpz_class> values = std::vector<mpz_class>();
            long tableOutputWidth = inputWidth;
            REPORT(DEBUG,"calculating Table " << i);

            for(long j = 0; j < pow(2, i); ++j) {
                long y = (j << (outputWidth - i)) + (1 << (outputWidth - i - 1));
                mpz_class inverse = calculateInverse(y);
                mpz_class lut_out;
                eval(lut_out, inverse);
                if(!monotoneIncreasing && lut_out.get_si() == y) {
                    ++inverse;
                }

                values.push_back(inverse);

                if(mpz_cmp(values.back().get_mpz_t(), inputWidth_ref.get_mpz_t()) >= 0) {
                    tableOutputWidth = inputWidth + 1;
                }
            }

            ComparatorTable *ct = new ComparatorTable(this, getTarget(), i, tableOutputWidth, values);
            addSubComponent(ct);

            string signal_table_out = declare(getTarget()->tableDelay(i, tableOutputWidth, true),
                                        join("ref", i), tableOutputWidth);

            string signal_comp_res = declare(getTarget()->adderDelay(tableOutputWidth) * i + getTarget()->tableDelay(i, tableOutputWidth, true) * (i-1),
                                             join("comp_res", i), 1);

            double table_in_delay = 0;
            if(i > 1) {
                table_in_delay = getTarget()->logicDelay(i);
            }

            string signal_table_in = declare(table_in_delay, join("table_input_", i), i);


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
            this->outPortMap(ct, "Y", signal_table_out);

            vhdl << this->instance(ct, join("ct", i));

            vhdl << tab << signal_comp_res << "(0) <= '1' when unsigned(i) " << comparison << " unsigned(" << signal_table_out << ") else '0';" << endl << endl;
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
    }


}//namespace