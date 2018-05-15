//
// Created by Viktor Schmidt on 2/9/18.
//

#include "MonotoneFunctionDiff.hpp"

using namespace std;
namespace flopoco {
    MonotoneFunctionDiff::MonotoneFunctionDiff(OperatorPtr parentOp, flopoco::Target *target, string functionString_, int inputWidth_, int outputWidth_)
            : FixMonotoneFunctionInterface(parentOp, target, functionString_, inputWidth_, outputWidth_) {

        srcFileName = "FixMonotoneFunctionDiff";

        ostringstream name;
        name << "FixMonotoneFunctionDiff" << inputWidth << "_" << outputWidth;
        setName(name.str());

        REPORT(INFO,"Declaration of FixMonotoneFunctionDiff \n");

        // more detailed message
        REPORT(DETAILED, "this operator has received two parameters " << inputWidth << " and " << outputWidth);

        build();
    }


    void MonotoneFunctionDiff::build() {
        mpz_class inputWidth_ref = mpz_class(1) << (inputWidth + 1);
        long tableOutputWidth = inputWidth + 1;

        string negate = monotoneIncreasing ? "not " : "";
//        string negate = "not ";
        vector<string> diff_signals(outputWidth);
        vector<string> ref_signals(outputWidth);
        vector<string> bit_signals(outputWidth);

        mpz_class inverse = mpz_class();
		vector<vector<mpz_class>> tables(outputWidth);
		vector<vector<mpz_class>> values(outputWidth);

        tables[0] = vector<mpz_class>();
        values[0] = vector<mpz_class>();

        inverse = calculateInverse(1 << (outputWidth - 1));
        mpz_class lut_out;

        eval(lut_out, inverse);

        if(!monotoneIncreasing && lut_out.get_si() == 1 << (outputWidth - 1)) {
            inverse += 1;
        }

        values[0].emplace_back(inverse);

        makeTwosComplement(inverse, tableOutputWidth);


        tables[0].emplace_back(inverse);

        bit_signals[0] = declare(monotoneIncreasing ? getTarget()->logicDelay(1) : 0, "output_0", 1);
        diff_signals[0] = declare(getTarget()->adderDelay(inputWidth + 1), "diff_0", tableOutputWidth);

        vhdl << tab << diff_signals[0]
             << " <= std_logic_vector(unsigned(i) + to_unsigned(2#" << inverse.get_str(2) << "#, "
             << tableOutputWidth << "));" << endl;

        vhdl << tab << bit_signals[0]
             << "(0) <= " << negate << diff_signals[0]
             << "(" << inputWidth<< ");" << endl << endl;

        // calculating table i
        for(int i = 1; i < outputWidth; ++i) {
            tables[i] = vector<mpz_class>();
            values[i] = vector<mpz_class>();
            REPORT(DEBUG,"calculating LUT " << i);

            for(int j = 0; j < pow(2, i); ++j) {
                long y = (j << (outputWidth - i)) + (1 << (outputWidth - i - 1));
                inverse = calculateInverse(y);

                if(!monotoneIncreasing && lut_out.get_si() == y) {
                    inverse += 1;
                }

                values[i].emplace_back(inverse);
                mpz_sub(inverse.get_mpz_t(), inverse.get_mpz_t(), values[i-1].at(j/2).get_mpz_t());

                makeTwosComplement(inverse, inputWidth + 1);

                tables[i].emplace_back(inverse);

                if(mpz_cmp(tables[i].back().get_mpz_t(), inputWidth_ref.get_mpz_t()) >= 0) {
                    tableOutputWidth = inputWidth + 2;
                }
            }

            ComparatorTable *ct = new ComparatorTable(this, getTarget(), i, tableOutputWidth, tables[i]);
            addSubComponent(ct);

            string signal_table_in = declare(join("table_input_", i), i);

            // creating the signals with delays
            double output_delay = monotoneIncreasing ? getTarget()->logicDelay(1) : 0;
            bit_signals[i] = declare(output_delay, join("output_", i), 1, false);

            double ref_delay = getTarget()->tableDelay(i, tableOutputWidth, true);
            ref_signals[i] = declare(ref_delay, join("ref_", i), tableOutputWidth);

            double diff_delay = getTarget()->adderDelay(tableOutputWidth);
            diff_signals[i] = declare(diff_delay, join("diff_", i), tableOutputWidth);

            // connecting table
            vhdl << tab << signal_table_in << " <= " << connectBits(bit_signals, 0, i-1) << ";" << endl << endl;
            this->inPortMap(ct, "X", signal_table_in);
            this->outPortMap(ct, "Y", ref_signals[i]);

            vhdl << this->instance(ct, join("ct", i)) << endl;
            vhdl << tab << diff_signals[i] << " <= std_logic_vector(unsigned(" << diff_signals[i-1] << ") + unsigned(" << ref_signals[i] << "));" << endl;
            vhdl << tab << bit_signals[i] << " <= " << negate << diff_signals[i] << "(" << tableOutputWidth-1 << ");" << endl;
        }

        vhdl << tab << "o <= " << connectBits(bit_signals, 0, outputWidth-1) << ";" << endl << endl;
    }

    OperatorPtr flopoco::MonotoneFunctionDiff::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
        string func;
        int inW, outW;
        UserInterface::parseString(args, "function", &func);
        UserInterface::parseInt(args, "inputWidth", &inW);
        UserInterface::parseInt(args, "outputWidth", &outW);
        return new MonotoneFunctionDiff(parentOp, target, func, inW, outW);
    }

    void MonotoneFunctionDiff::registerFactory() {
        UserInterface::add("MonotoneFunctionDiff", // name
                           "Generates a function.", // description, string
                           "Miscellaneous", // category, from the list defined in UserInterface.cpp
                           "",
                           "function(string)=x: Algorithm Type: normal, diff, lut;\
                            inputWidth(int)=16: Input bit count; \
                            outputWidth(int)=8: Output bit count",
                           "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developer manual in the doc/ directory of FloPoCo.",
                           MonotoneFunctionDiff::parseArguments
        );
    }
}
