#ifndef CONVOLUTIONALCORESIMPLE_H
#define CONVOLUTIONALCORESIMPLE_H
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
	class ConvolutionalCoreSimple : public Operator {

    public:
        ConvolutionalCoreSimple(Target* target, unsigned int wordSize_, unsigned int fraction_, unsigned int weightWordSize_, unsigned int weightFraction_, unsigned int size_, vector <double> weights_, bool useAdderTree_=true, string roundingType_="Truncation", string id_="0");

		// destructor
        ~ConvolutionalCoreSimple() {}
        unsigned int getOutputWordSize();
        unsigned int getOutputFraction();

    private:
        unsigned int outputWordSize;
        unsigned int outputFraction;
        unsigned int wordSize;
        unsigned int fraction;
        unsigned int weightWordSize;
        unsigned int weightFraction;
        unsigned int size;
		vector <double> weights;
        bool useAdderTree;
		string id;
        vector <string> inputNames;

        string roundingType;
        void roundOutput(unsigned int wordSizeFrom, unsigned int fractionFrom, string round="Truncation");
	};


}//namespace

#endif
