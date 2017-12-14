#ifndef PADDINGGENERATOR_H
#define PADDINGGENERATOR_H
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
	class PaddingGenerator : public Operator {

    public:
        PaddingGenerator(Target* target, unsigned int wordSize_, unsigned int windowSize_, unsigned int horizontalSize_, unsigned int verticalSize_, int padTop_=1, unsigned int stride_=1, string padType_="Zero", int padBot_=-1, int padLeft_=-1, int padRight_=-1, bool genValidFinished_=true);

		// destructor
        ~PaddingGenerator() {}

        string getPaddingValue(string name);

        vector <vector <string>> inputNames;
        vector <vector <string>> outputNames;

	private:
		unsigned int wordSize;
        string padType;
        unsigned int windowSize;
        unsigned int horizontalSize;
        unsigned int verticalSize;
        unsigned int stride;
        int padTop;
        int padBot;
        int padLeft;
        int padRight;
        bool genValidFinished;

        unsigned int numberOfInputs;
	};


}//namespace

#endif
