#ifndef WINDOWSHIFTREGISTER
#define WINDOWSHIFTREGISTER
/* Each Operator declared within the flopoco framework has 
   to inherit the class Operator and overload some functions listed below*/
#include "Operator.hpp"

/* This file contains a lot of useful functions to manipulate vhdl */
#include "utils.hpp"
#include <string>
#include <vector>

/*  All flopoco operators and utility functions are declared within
    the flopoco namespace.
    You have to use flopoco:: or using namespace flopoco in order to access these
    functions.
*/

namespace flopoco {

	// new operator class declaration
    class WindowShiftRegister : public Operator {

    public:
        WindowShiftRegister(Target* target, unsigned int wordSize_, unsigned int windowSize_, unsigned int rowLength_);

		// destructor
        ~WindowShiftRegister() {}
        vector <string> outputNames;

	private:
        unsigned int numberOfFlipFlops;
        unsigned int wordSize;
        unsigned int windowSize;
        unsigned int rowLength;
        unsigned int stride;
	};


}//namespace

#endif
