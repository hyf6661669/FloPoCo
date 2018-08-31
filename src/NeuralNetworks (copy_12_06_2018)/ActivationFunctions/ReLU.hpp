#ifndef RELU_H
#define RELU_H
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
	class ReLU : public Operator {

    public:
        ReLU(Target* target, unsigned int wordSize_);

		// destructor
        ~ReLU() {}

	private:
		unsigned int wordSize;
	};


}//namespace

#endif
