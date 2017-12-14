#ifndef POOLINGLAYER_H
#define POOLINGLAYER_H
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
    class PoolingLayer : public Operator {

    public:
        PoolingLayer(Target* target, unsigned int wordSize_, unsigned int horizontalSize_, unsigned int verticalSize_, unsigned int numberOfFeatures_, bool parallel_, int paddingTop_, unsigned int windowSize_=2, unsigned int stride_=1, int paddingBot_=-1, int paddingLeft_=-1, int paddingRight_=-1, string paddingType_="Zero");

		// destructor
        ~PoolingLayer() {}

        unsigned int getNumberOfInstances();

    private:
	unsigned int wordSize;
        unsigned int horizontalSize;
        unsigned int verticalSize;
        unsigned int numberOfFeatures;
        bool parallel;
        int paddingTop=-1;
        unsigned int windowSize;
        unsigned int stride;
        int paddingBot;
        int paddingLeft;
        int paddingRight;
        string paddingType;
	};


}//namespace

#endif
