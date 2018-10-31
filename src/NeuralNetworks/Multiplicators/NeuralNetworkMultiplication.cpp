//
// Created by nfiege on 28/06/18.
//
// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "NeuralNetworkMultiplication.hpp"

#include "SimpleMultiplication.cpp"

namespace flopoco {
    NeuralNetworkMultiplication::NeuralNetworkMultiplication(Target *target, int wIn1_, int wIn2_, multiplicationType mType) :
            Operator(target), wIn1(wIn1_), wIn2(wIn2_) {
        useNumericStd();
        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName = "NeuralNetworkMultiplication";

        // author
        setCopyrightString("Nicolai Fiege, 2018");

        // definition of the name of the operator
        ostringstream name;
        name << "NeuralNetworkMultiplication_" << wIn1 << "_" << wIn2 << "_" << mType;
        setName(name.str());

        this->wOut = wIn1 + wIn2;


        addInput("X", wIn1);
        addInput("Y", wIn2);

        switch (mType)
        {
            case simple: this->buildSimpleMultiplier(target); break;
            default: THROWERROR("Undefined multiplication type!"); break;
        }

        this->syncCycleFromSignal("R_temp");
        vhdl << "R <= R_temp;" << endl;

        addOutput("R", wOut);
    }

    void NeuralNetworkMultiplication::buildSimpleMultiplier(Target* target)
    {
        SimpleMultiplication* op = new SimpleMultiplication(target,this->wIn1,this->wIn2);
        addSubComponent(op);
        inPortMap(op,"X","X");
        inPortMap(op,"Y","Y");
        outPortMap(op,"R","R_temp",true);
        vhdl << instance(op,"Simple_Multiplication_instance");
    }
}