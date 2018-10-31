//
// Created by nfiege on 28/06/18.
//
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

#include "SimpleMultiplication.hpp"


namespace flopoco{
    SimpleMultiplication::SimpleMultiplication(Target *target, int wIn1, int wIn2) : Operator(target)
    {
        useNumericStd();
        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName = "SimpleMultiplication";

        // author
        setCopyrightString("Nicolai Fiege, 2018");

        // definition of the name of the operator
        ostringstream name;
        name << "SimpleMultiplication_" << wIn1 << "_" << wIn2;
        setName(name.str());


        addInput("X", wIn1);
        addInput("Y", wIn2);
        this->nextCycle(); // input register

        addOutput("R", wIn1+wIn2);

        vhdl << declare("R_temp",wIn1+wIn2) << " <= std_logic_vector(signed(X) * signed(Y));" << endl;
        this->nextCycle(); // pipelined multiplication (dirty version...)

        vhdl << "R <= R_temp;" << endl;
    }
}