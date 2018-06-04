#ifndef MAX2INPUT_H
#define MAX2INPUT_H
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
	class Max2Input : public Operator {

    public:
        Max2Input(Target* target, unsigned int wordSize_, bool sign_);

		// destructor
        ~Max2Input() {}

	private:
		unsigned int wordSize;
        bool sign;
	};


}//namespace

#endif
