//
// Created by Viktor Schmidt on 3/5/18.
//

#ifndef FLOPOCO_FIXMONOTONEFUNCTIONINTERFACE_HPP
#define FLOPOCO_FIXMONOTONEFUNCTIONINTERFACE_HPP

#include "Operator.hpp"
#include "utils.hpp"
#include <string>
#include <iomanip>
#include <sollya.h>
#include <sstream>
#include "gmp.h"
#include "mpfr.h"

namespace flopoco {
    class FixMonotoneFunctionInterface : public Operator {
    protected:
        int inputWidth;
        int outputWidth;

        bool monotoneIncreasing;

        string functionString;
        sollya_obj_t fS;


        void eval(mpz_class &r, mpz_class x) const;
        bool checkMonotoneIncreasing();
        void makeTwosComplement(mpz_class &r, int bitCount);
        virtual void build() = 0;

    public:
        FixMonotoneFunctionInterface(OperatorPtr parentOp, Target *target, string functionString_, int inputWidth_, int outputWidth_);
        //~FixMonotoneFunctionInterface() {};

        void emulate(TestCase *tc) override;
        void buildStandardTestCases(TestCaseList *tcl) override;

    };
}

#endif //FLOPOCO_FIXMONOTONEFUNCTIONINTERFACE_HPP
