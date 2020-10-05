#include <iostream>
#include <sstream>
#include <string>
#include "gmp.h"
#include "mpfr.h"
#include <vector>
#include <gmpxx.h>
#include <stdio.h>
#include <stdlib.h>
#include "ModReduction/PseudoCompressor.hpp"


using namespace std;
namespace flopoco {

    PseudoCompressor::PseudoCompressor(Operator *parentOp, Target *target, int weight, int modulus) : Compressor(
            parentOp, target) {
        setCopyrightString("Andreas Boettcher");

        this->_weight = weight;
        this->_modulus = modulus;

        ostringstream name;
        name << "Pseudo_Compressor_Weight_" << _weight << "_for_Modulus_" << _modulus;
        setNameWithFreqAndUID(name.str());

        heights.resize(_weight, 0);
        heights[_weight - 1] = 1;               //the pseudocompressor takes 1 bit from the column with the provided weight

/*        for (unsigned i = 0; i < heights.size(); i++) {
            if (heights[heights.size() - i - 1] > 0)
                addInput(join("X", i), heights[heights.size() - i - 1]);
        }*/
        REPORT(DEBUG, "add pseudocompressor input with weight " << _weight << " for modulus of " << modulus);
        addInput(join("X", _weight-1), 1);

        unsigned long long place_value = 1;
        for (int j = 1; j <= _weight; j++)  place_value *= 2;
        REPORT(DEBUG, "place value is " << place_value );


        unsigned long long offline_mod = place_value % modulus;
        REPORT(DEBUG, "modulus contribution " << offline_mod );
        for(unsigned i = 1, j = 0; i < offline_mod; i *= 2, j++){       //as long as there are more bits 1 in the modulus
            cout << " i is " << i << " j is" << j << endl;
            if(offline_mod & i){
                addOutput(join("R", j), 1);
                REPORT(DEBUG, "add pseudocompressor output with weight " << i << " for modulus of " << modulus);
                vhdl << tab << join("R", j) << " <= " << join("X", _weight-1) << ";" << endl;
            }

        }

    }

    PseudoCompressor::~PseudoCompressor() {
    }

    OperatorPtr PseudoCompressor::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
        int weight, modulus;
        UserInterface::parseInt(args, "weight", &weight);
        UserInterface::parseInt(args, "mod", &modulus);

        return new PseudoCompressor(parentOp, target, weight, modulus);
    }

    void PseudoCompressor::registerFactory() {
        UserInterface::add("PseudoCompressor", // name
                           "1:mod pseudo compressor used for modulo reduction.", // description, string
                           "Primitives", // category, from the list defined in UserInterface.cpp
                           "",
                           "weight(int): n of 2^n represents the place value the modulo is calculated offline for; \
                        mod(int)=: Specifies the modulus for the calculation",
                           "",
                           PseudoCompressor::parseArguments
        );
    }

}
