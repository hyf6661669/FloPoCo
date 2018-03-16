//
// Created by Viktor Schmidt.
//

#include "MonotoneFunction.hpp"

using namespace std;
namespace flopoco {
    MonotoneFunction::MonotoneFunction(OperatorPtr parentOp, Target* target, int inputWidth_, int outputWidth_)
            : FixMonotoneFunctionInterface(parentOp, target, inputWidth_, outputWidth_) {
        srcFileName="FixMonotoneFunction";
        //useNumericStd();

        // definition of the name of the operator
        ostringstream name;
        name << "FixMonotoneFunction" << inputWidth << "_" << outputWidth;
        setName(name.str());

        REPORT(INFO,"Declaration of FixMonotoneFunction \n");

        // more detailed message
        REPORT(DETAILED, "this operator has received two parameters " << inputWidth << " and " << outputWidth);

        // debug message for developper
        REPORT(DEBUG,"debug of FixMonotoneFunction");
        //sollya_lib_set_prec(sollya_lib_constant_from_int(65536));
//        sollya_lib_parse_string("[0;1]");
//        fS= sollya_lib_parse_string("1 - sin(3.14/4*x) - 1/256;");

        build();
    };


    mpz_class MonotoneFunction::calculateInverse(int y) {
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

        return lut_in;
    }


    OperatorPtr MonotoneFunction::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
        int param0, param1;
        UserInterface::parseInt(args, "inputWidth", &param0);
        UserInterface::parseInt(args, "outputWidth", &param1);
        return new MonotoneFunction(parentOp, target, param0, param1);
    }

    void MonotoneFunction::registerFactory(){
        UserInterface::add("MonotoneFunction", // name
                           "Generates a function.", // description, string
                           "Miscellaneous", // category, from the list defined in UserInterface.cpp
                           "", //seeAlso
                // Now comes the parameter description string.
                // Respect its syntax because it will be used to generate the parser and the docs
                // Syntax is: a semicolon-separated list of parameterDescription;
                // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
                           "inputWidth(int)=16: Input bit count; \
                        outputWidth(int)=8: Output bit count",
                // More documentation for the HTML pages. If you want to link to your blog, it is here.
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