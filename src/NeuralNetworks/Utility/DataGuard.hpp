#ifndef DATAGUARD_H
#define DATAGUARD_H
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
	class DataGuard : public Operator {

    public:
        DataGuard(Target* target, int wordSize_, int wordSizeGuard_, unsigned int numberToCompare_);

		// destructor
        ~DataGuard() {}

	private:
		unsigned int wordSize;
		int wordSizeGuard;
		unsigned int numberToCompare;
	};


}//namespace

#endif // DATAGUARD_H
