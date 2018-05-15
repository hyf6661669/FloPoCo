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
        TestCase *tc;
        mpz_class max_input = (mpz_class(1) << inputWidth) - 1;

        tc = new TestCase(this);
        tc->addInput("i", 0);
        emulate(tc);
        tcl->add(tc);

        for(long y = 0; y < pow(2, outputWidth); ++y) {
            mpz_class inverse = calculateInverse(y);

            mpz_class actual_y;
            eval(actual_y, inverse);

            // only add test case if the actual output matches the input
            // (not always the case for functions that don't use the full output range)
            // also don't add 0 or the last value because they are added separately
            if(actual_y.get_si() == y && inverse.get_si() != 0 && inverse != max_input) {
                tc = new TestCase(this);
                tc->addInput("i", inverse);
                emulate(tc);
                tcl->add(tc);

                if(monotoneIncreasing) {
                    --inverse;
                }
                else {
                    ++inverse;
                }

                tc = new TestCase(this);
                tc->addInput("i", inverse);
                emulate(tc);
                tcl->add(tc);
            }
        }

        tc = new TestCase(this);
        tc->addInput("i", max_input);
        emulate(tc);
        tcl->add(tc);
    }

    mpz_class FixMonotoneFunctionInterface::calculateInverse(long y) {
        REPORT(FULL,"calculateInverse looking for x at f(x)=" << y);

        mpz_class inverse = mpz_class(1) << (inputWidth - 1);
        mpz_class lut_out;
        long sign = monotoneIncreasing ? 1 : -1;

        for(int i = inputWidth - 2; i >= 0; --i) {
            eval(lut_out, inverse);

            if(lut_out.get_si() >= y) {
                inverse -= sign * mpz_class(1) << i;
            }
            else {
                inverse += sign * mpz_class(1) << i;
            }
        }

        eval(lut_out, inverse);

        if(lut_out.get_si() >= y) {
            inverse -= sign;
        }

        eval(lut_out, inverse);

        if(lut_out.get_si() < y) {
            inverse += sign;
        }

        if(inverse > (mpz_class(1) << inputWidth)) {
            inverse = mpz_class(1) << inputWidth;
            REPORT(INFO, "Clipping occured.");
        }

        if(inverse < 0) {
            inverse = mpz_class(0);
            REPORT(INFO, "Clipping occured.");
        }

        return inverse;
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

    string FixMonotoneFunctionInterface::connectBits(vector<string> signals, int start, int end) {
        ostringstream connect;

        for(int i = start; i <= end; ++i) {
            connect << signals.at(i);

            if(i < end) {
                connect << " & ";
            }
        }

        return connect.str();
    }

}