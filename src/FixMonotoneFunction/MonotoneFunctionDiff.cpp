//
// Created by Viktor Schmidt on 2/9/18.
//

#include "MonotoneFunctionDiff.hpp"

using namespace std;
namespace flopoco {
    MonotoneFunctionDiff::MonotoneFunctionDiff(flopoco::Target *target, int inputWidth_, int outputWidth_)
            : FixMonotoneFunctionInterface(target, inputWidth_, outputWidth_) {

        srcFileName = "MonotoneFunctionDiff";
        //useNumericStd();

        // definition of the name of the operator
        ostringstream name;
        name << "MonotoneFunctionDiff" << inputWidth << "_" << outputWidth;
        setName(name.str());

        //setCopyrightString("Viktor Schmidt 2018");

        // declaring inputs
        //addInput("i" , inputWidth);
//        addFullComment(" addFullComment for a large comment ");
//        addComment("addComment for small left-aligned comment");

        // declaring output
        //addOutput("o" , outputWidth);
        // basic message
        REPORT(INFO,"Declaration of MonotoneFunctionDiff \n");

        // more detailed message
        REPORT(DETAILED, "this operator has received two parameters " << inputWidth << " and " << outputWidth);

        // debug message for developper
        REPORT(DEBUG,"debug of MonotoneFunctionDiff");
        //sollya_lib_set_prec(sollya_lib_constant_from_int(65536));
//        sollya_lib_parse_string("[0;1]");
//        fS= sollya_lib_parse_string("1 - sin(3.14/4*x) - 1/256;");
//
//        mpz_class in_0, out_0, in_max, out_max;
//        in_0 = mpz_class(0);
//        in_max = pow(2, inputWidth) - 1;
//        eval(out_0, in_0);
//        eval(out_max, in_max);
//
//        monotoneIncreasing = mpz_cmp(out_max.get_mpz_t(), out_0.get_mpz_t()) >= 0;

        build();
    }
//    void MonotoneFunctionDiff::emulate(TestCase * tc) {
//        mpz_class sx = tc->getInputValue("i");
//
//        mpz_class sr;
//        eval(sr, sx);
//
//        tc->addExpectedOutput("o",sr);
//    }


//    void MonotoneFunctionDiff::buildStandardTestCases(TestCaseList * tcl) {
//        // please fill me with regression tests or corner case tests!
//    }


//    void MonotoneFunctionDiff::eval(mpz_class& r, mpz_class x) const
//    {
//        mpfr_t mpX, mpR;
//        mpfr_init2(mpX, inputWidth+2);
//        mpfr_init2(mpR, outputWidth*3);
//        sollya_lib_set_prec(sollya_lib_constant_from_int(inputWidth*2));
//
//        // scaling to [0, 1]
//        mpfr_set_z(mpX, x.get_mpz_t(), GMP_RNDN);
//        mpfr_div_2si(mpX, mpX, inputWidth, GMP_RNDN);
//
//        sollya_lib_evaluate_function_at_point(mpR, fS, mpX, NULL);
//
//        mpfr_mul_2si(mpR, mpR, outputWidth, GMP_RNDN);
//        mpfr_get_z(r.get_mpz_t(), mpR, GMP_RNDN);
//
//        REPORT(FULL,"f("<< mpfr_get_d(mpX, GMP_RNDN) << ") = " <<mpfr_get_d(mpR, GMP_RNDN));
//
//        mpfr_clear(mpX);
//        mpfr_clear(mpR);
//    }

    mpz_class MonotoneFunctionDiff::calculateInverse(int y) {
        REPORT(DEBUG,"calculateInverse looking for x at f(x)=" << y);
        mpz_class y_real = y;
        int ref = (int)pow(2, inputWidth) - 1;
        int nextOffset = (int)pow(2, inputWidth - 1);
        mpz_class lut_in(ref), lut_out;
        int sign = monotoneIncreasing ? 1 : -1;

        for(int i = 0; i <= inputWidth; i++) {
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

        if(mpz_cmp(lut_out.get_mpz_t(), y_real.get_mpz_t()) == sign) {
            ref += sign;
        }

        lut_in = mpz_class(ref);

        return lut_in;
    }


    void MonotoneFunctionDiff::build() {
        string negate = monotoneIncreasing ? "not " : "";
        declare("output", outputWidth);
        mpz_class r = mpz_class();
        vector<mpz_class> tables[outputWidth];
        vector<mpz_class> values[outputWidth];

        tables[0] = vector<mpz_class>();
        values[0] = vector<mpz_class>();

        r = calculateInverse(1 << (outputWidth - 1));
        values[0].push_back(mpz_class(r));

        if(mpz_cmp_d(r.get_mpz_t(), 0) > 0) {
            mpz_abs(r.get_mpz_t(), r.get_mpz_t());
            r -= mpz_class(1) << (inputWidth + 1);
            mpz_abs(r.get_mpz_t(), r.get_mpz_t());
        }
        else {
            mpz_abs(r.get_mpz_t(), r.get_mpz_t());
        }

        tables[0].push_back(mpz_class(r));

        string diffSignals[outputWidth];
        diffSignals[0] = declare(join("diff", 0), inputWidth + 1);

        vhdl << tab << diffSignals[0] << " <= std_logic_vector(unsigned(i) + to_unsigned(2#" << r.get_str(2) << "#, " << inputWidth + 1 << "));" << endl;
        vhdl << tab << "output(" << outputWidth - 1 << ") <= " << negate << diffSignals[0] << "(" << inputWidth<< ");" << endl << endl;

        for(int i = 1; i < outputWidth; ++i) {
            tables[i] = vector<mpz_class>();
            values[i] = vector<mpz_class>();
            REPORT(DEBUG,"calculating LUT " << i);
            for(int j = 0; j < pow(2, i); ++j) {
                //r = mpz_class(0);
                //REPORT(DEBUG,r.get_str(10) << " == " << r.get_str(2));
                int v = (j << (outputWidth - i)) + (1 << (outputWidth - i - 1));
                r = calculateInverse(v);
                values[i].push_back(mpz_class(r));
                //mpz_add(r.get_mpz_t(), calculateInverse(v).get_mpz_t(), tables[i-1].at(j/2).get_mpz_t());
                mpz_sub(r.get_mpz_t(), r.get_mpz_t(), values[i-1].at(j/2).get_mpz_t());
                REPORT(DEBUG,"sub result: " << values[i].at(j).get_str(2) << " - " << values[i-1].at(j/2).get_str(2) << " = " << r.get_str(2));


                if(mpz_cmp_d(r.get_mpz_t(), 0) > 0) {
                    //REPORT(DEBUG,r.get_str(10) << " = " << r.get_str(2));
                    mpz_abs(r.get_mpz_t(), r.get_mpz_t());
                    //mpz_com(r.get_mpz_t(), r.get_mpz_t());
                    r -= mpz_class(1) << (inputWidth + 1);
                    mpz_abs(r.get_mpz_t(), r.get_mpz_t());

                    //REPORT(DEBUG,r.get_str(10) << " = " << r.get_str(2));
                    //mpz_neg(r.get_mpz_t(), r.get_mpz_t());
                }
                else {
                    mpz_abs(r.get_mpz_t(), r.get_mpz_t());
                }

                tables[i].push_back(mpz_class(r));
            }

            ComparatorTable *ct = new ComparatorTable(this, getTarget(), i, inputWidth + 1, tables[i]);
            addSubComponent(ct);

            string signal = declare(join("ref", i), inputWidth + 1);
            diffSignals[i] = declare(join("diff", i), inputWidth + 1);
            vhdl << tab << "ct" << i << " : " << ct->getName() << " port map(X => output(" << outputWidth - 1 << " downto " << outputWidth - i << "), Y => " << signal << ", clk => clk, rst => rst);" << endl;
            vhdl << tab << diffSignals[i] << " <= std_logic_vector(unsigned(" << diffSignals[i-1] << ") + unsigned(" << signal << "));" << endl;
            vhdl << tab << "output(" << outputWidth - i - 1 << ") <= " << negate << diffSignals[i] << "(" << inputWidth << ");" << endl << endl;
        }

        vhdl << tab << "o <= output;" << endl;
    }

    OperatorPtr flopoco::MonotoneFunctionDiff::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
        int param0, param1;
        UserInterface::parseInt(args, "inputWidth", &param0);
        UserInterface::parseInt(args, "outputWidth", &param1);
        return new MonotoneFunctionDiff(target, param0, param1);
    }

    void MonotoneFunctionDiff::registerFactory() {
        UserInterface::add("MonotoneFunctionDiff", // name
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
                           MonotoneFunctionDiff::parseArguments
        );
    }
}
