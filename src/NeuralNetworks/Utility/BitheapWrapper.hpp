//
// Created by nfiege on 23/05/18.
//


#ifndef FLOPOCO_BITHEAPWRAPPER_HPP
#define FLOPOCO_BITHEAPWRAPPER_HPP
/* Each Operator declared within the flopoco framework has
   to inherit the class Operator and overload some functions listed below*/
#include "Operator.hpp"

/* This file contains a lot of useful functions to manipulate vhdl */
#include "utils.hpp"
#include <string>

/*  All flopoco operators and utility functions are declared within
    the flopoco namespace.
    You have to use flopoco:: or using namespace flopoco in order to access these
    functions.
*/

namespace flopoco {

    // new operator class declaration
    class BitheapWrapper : public Operator {

    public:
        BitheapWrapper(Target* target, vector<int> wIn, int numberOfInputs, vector<bool> signs, vector<int> weights);

        // destructor
        ~BitheapWrapper() {}
        unsigned int getOutputWordSize() const {return this->wOut;}

    private:
        unsigned int wOut;
    };


}//namespace

#endif // FLOPOCO_BITHEAPWRAPPER_HPP