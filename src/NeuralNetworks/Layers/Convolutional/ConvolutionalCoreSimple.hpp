#ifndef CONVOLUTIONALCORESIMPLE_H
#define CONVOLUTIONALCORESIMPLE_H
/* Each Operator declared within the flopoco framework has 
   to inherit the class Operator and overload some functions listed below*/
#include "Operator.hpp"

/* This file contains a lot of useful functions to manipulate vhdl */
#include "utils.hpp"
#include <string>
#include <vector>


namespace flopoco {

	// new operator class declaration
	class ConvolutionalCoreSimple : public Operator {

    public:
        ConvolutionalCoreSimple(Target* target, unsigned int wordSize_, unsigned int fraction_, unsigned int weightWordSize_, unsigned int weightFraction_, unsigned int size_, vector <double> weights_, bool useBitHeap_=true, char roundingType_=0x00, string id_="0");

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
        bool useBitHeap;
		string id;
        vector <string> inputNames;

        char roundingType; // 0x00: Truncation, 0x01: Saturation
        void roundOutput(Target* target, unsigned int wordSizeFrom, unsigned int fractionFrom, char round = 0x00);
	};


}//namespace

#endif
