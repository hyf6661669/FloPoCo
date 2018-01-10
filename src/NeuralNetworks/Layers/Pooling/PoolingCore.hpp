#ifndef POOLINGCORE_H
#define POOLINGCORE_H
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
	class PoolingCore : public Operator {

    public:
        PoolingCore(Target* target, unsigned int wordSize_, unsigned int windowSize_);

		// destructor
        ~PoolingCore() {}

	private:
        unsigned int wordSize;
        unsigned int windowSize;
		vector <string> inputNames;
	};


}//namespace

#endif
