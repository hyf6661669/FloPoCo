//
// Created by Viktor Schmidt on 2/9/18.
//

#include "MonotoneFunctionDiff.hpp"

using namespace std;
namespace flopoco {
    MonotoneFunctionDiff::MonotoneFunctionDiff(OperatorPtr parentOp, flopoco::Target *target, string functionString_, int inputWidth_, int outputWidth_)
            : FixMonotoneFunctionInterface(parentOp, target, functionString_, inputWidth_, outputWidth_) {

        srcFileName = "MonotoneFunctionDiff";

        ostringstream name;
        name << "MonotoneFunctionDiff" << inputWidth << "_" << outputWidth;
        setName(name.str());

        REPORT(INFO,"Declaration of MonotoneFunctionDiff \n");

        // more detailed message
        REPORT(DETAILED, "this operator has received two parameters " << inputWidth << " and " << outputWidth);

        build();
    }


    mpz_class MonotoneFunctionDiff::calculateInverse(int y) {
        REPORT(DEBUG,"calculateInverse looking for x at f(x)=" << y);
        mpz_class y_real = y;
        int ref = (int)pow(2, inputWidth) - 1;
        int nextOffset = (int)pow(2, inputWidth - 1);
        mpz_class lut_in(ref), lut_out;
        int sign = monotoneIncreasing ? 1 : -1;

        for(int i = 0; i < inputWidth; i++) {
            lut_in = mpz_class(ref);
            eval(lut_out, lut_in);

            if(mpz_cmp(lut_out.get_mpz_t(), y_real.get_mpz_t()) >= 0) {
                ref -= sign * nextOffset;
            }
            else {
                ref += sign * nextOffset;
            }

            nextOffset /= 2;
        }

        lut_in = mpz_class(ref);
        eval(lut_out, lut_in);

        if(mpz_cmp(lut_out.get_mpz_t(), y_real.get_mpz_t()) < 0) {
                ref += sign;
        }

        lut_in = mpz_class(ref);

        if(lut_in.get_si() > pow(2, inputWidth) - 1) {
            REPORT(DEBUG,"inverse too big: f(" << lut_in.get_str(2) << ")=" << y);
        }

        return lut_in;
    }


    void MonotoneFunctionDiff::build() {
        string negate = monotoneIncreasing ? "not " : "";
        double delay = 0;

        for(int x = 0; x < outputWidth; ++x) {
            delay += getTarget()->adderDelay(inputWidth+1);
            delay += getTarget()->tableDelay(x, inputWidth+1, true);
        }

        declare(delay, "output", outputWidth);
        mpz_class r = mpz_class();
		vector<vector<mpz_class>> tables(outputWidth);
		vector<vector<mpz_class>> values(outputWidth);

        tables[0] = vector<mpz_class>();
        values[0] = vector<mpz_class>();

        r = calculateInverse(1 << (outputWidth - 1));

        if(!monotoneIncreasing) {
            r += 1;
        }

        values[0].emplace_back(r);

        makeTwosComplement(r, inputWidth + 1);


        tables[0].emplace_back(r);

		vector<string> diffSignals(outputWidth);
        diffSignals[0] = declare(getTarget()->adderDelay(inputWidth + 1), join("diff", 0), inputWidth + 1);

        vhdl << tab << diffSignals[0] << " <= std_logic_vector(unsigned(i) + to_unsigned(2#" << r.get_str(2) << "#, " << inputWidth + 1 << "));" << endl;
        vhdl << tab << "output(" << outputWidth - 1 << ") <= " << negate << diffSignals[0] << "(" << inputWidth<< ");" << endl << endl;

        for(int i = 1; i < outputWidth; ++i) {
            tables[i] = vector<mpz_class>();
            values[i] = vector<mpz_class>();
            REPORT(DEBUG,"calculating LUT " << i);

            for(int j = 0; j < pow(2, i); ++j) {
                int v = (j << (outputWidth - i)) + (1 << (outputWidth - i - 1));
                r = calculateInverse(v);

                if(!monotoneIncreasing) {
                    r += 1;
                }

                values[i].emplace_back(r);
                mpz_sub(r.get_mpz_t(), r.get_mpz_t(), values[i-1].at(j/2).get_mpz_t());
                REPORT(DEBUG,"sub result: " << values[i].at(j).get_str(2) << " - " << values[i-1].at(j/2).get_str(2) << " = " << r.get_str(2));

                makeTwosComplement(r, inputWidth + 1);

                tables[i].emplace_back(r);
            }

            ComparatorTable *ct = new ComparatorTable(this, getTarget(), i, inputWidth + 1, tables[i]);
            addSubComponent(ct);

            string signal = declare(getTarget()->tableDelay(i, inputWidth + 1, true), join("ref", i), inputWidth + 1);
            diffSignals[i] = declare(getTarget()->adderDelay(inputWidth + 1), join("diff", i), inputWidth + 1);
            //string signal = declare(join("ref", i), inputWidth + 1);
            //diffSignals[i] = declare(join("diff", i), inputWidth + 1);

            ostringstream outputPort;
            outputPort << "output(" << outputWidth - 1 << " downto " << outputWidth - i << ")";
            this->inPortMap(ct, "X", outputPort.str());
            this->outPortMap(ct, "Y", signal);



            vhdl << this->instance(ct, join("ct", i));
            vhdl << tab << diffSignals[i] << " <= std_logic_vector(unsigned(" << diffSignals[i-1] << ") + unsigned(" << signal << "));" << endl;
            vhdl << tab << "output(" << outputWidth - i - 1 << ") <= " << negate << diffSignals[i] << "(" << inputWidth << ");" << endl << endl;
        }

        vhdl << tab << "o <= output;" << endl;
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
