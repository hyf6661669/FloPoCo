//
// Created by Viktor Schmidt on 3/5/18.
//

#include "FixMonotoneFunctionInterface.hpp"

using namespace std;
namespace flopoco {
    FixMonotoneFunctionInterface::FixMonotoneFunctionInterface(OperatorPtr parentOp, Target *target, string functionString_, int inputWidth_, int outputWidth_) :
            Operator(parentOp, target), functionString(functionString_), inputWidth(inputWidth_), outputWidth(outputWidth_) {

        useNumericStd();

        setCopyrightString("Viktor Schmidt 2018");

        addInput("i" , inputWidth);
        addOutput("o" , outputWidth);

        sollya_lib_parse_string("[0;1]");
        //fS= sollya_lib_parse_string("sin(3.14/2*x) * 255/256;");
        //fS = sollya_lib_parse_string("1 - sin(3.14/4*x) - 1/256;");
        fS = sollya_lib_parse_string(functionString.c_str());

        monotoneIncreasing = checkMonotoneIncreasing();

        if(monotoneIncreasing) {
            REPORT(DEBUG, "Function is monotone increasing.");
            vhdl << endl << "-- Implementing monotone increasing function f(x) = " << functionString << endl << endl;
        }
        else {
            REPORT(DEBUG, "Function is monotone decreasing.");
            vhdl << endl << "-- Implementing monotone decreasing function f(x) = " << functionString << endl << endl;
        }
    }

    void FixMonotoneFunctionInterface::emulate(TestCase *tc) {
        mpz_class sx = tc->getInputValue("i");

        mpz_class sr;
        eval(sr, sx);

        tc->addExpectedOutput("o", sr);
    }


    void FixMonotoneFunctionInterface::buildStandardTestCases(TestCaseList *tcl) {
        // please fill me with regression tests or corner case tests!
    }


    void FixMonotoneFunctionInterface::eval(mpz_class &r, mpz_class x) const {
        mpfr_t mpX, mpR;
        mpfr_init2(mpX, inputWidth * 2);
        mpfr_init2(mpR, outputWidth * 3);
        sollya_lib_set_prec(sollya_lib_constant_from_int(inputWidth * 2));

        // scaling to [0, 1]
        mpfr_set_z(mpX, x.get_mpz_t(), GMP_RNDN);
        mpfr_div_2si(mpX, mpX, inputWidth, GMP_RNDN);

        sollya_lib_evaluate_function_at_point(mpR, fS, mpX, NULL);

        mpfr_mul_2si(mpR, mpR, outputWidth, GMP_RNDN);
        mpfr_get_z(r.get_mpz_t(), mpR, GMP_RNDN);

        //REPORT(FULL, "f(" << mpfr_get_d(mpX, GMP_RNDN) << ") = " << mpfr_get_d(mpR, GMP_RNDN));

        mpfr_clear(mpX);
        mpfr_clear(mpR);
    }

    bool FixMonotoneFunctionInterface::checkMonotoneIncreasing() {
        mpz_class in_0, out_0, in_max, out_max;

        in_0 = mpz_class(0);
        in_max = pow(2, inputWidth) - 1;

        eval(out_0, in_0);
        eval(out_max, in_max);

        return mpz_cmp(out_max.get_mpz_t(), out_0.get_mpz_t()) >= 0;
    }

    void FixMonotoneFunctionInterface::makeTwosComplement(mpz_class &r, int bitCount) {
        if(mpz_cmp_d(r.get_mpz_t(), 0) > 0) {
            mpz_abs(r.get_mpz_t(), r.get_mpz_t());
            r -= mpz_class(1) << bitCount;
            mpz_abs(r.get_mpz_t(), r.get_mpz_t());
        }
        else {
            mpz_abs(r.get_mpz_t(), r.get_mpz_t());
        }
    }

}