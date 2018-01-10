#ifndef OUTPUTLAYER_H
#define OUTPUTLAYER_H
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
#include "NeuralNetworks/Layers/Layer.hpp"

namespace flopoco {

	// new operator class declaration
    class OutputLayer : public Layer {

    public:
        OutputLayer(Target* target, int howMany_, int wordSize_);

		// destructor
        ~OutputLayer() {}

        virtual string getOutputSignalName(int feature) override;
        virtual string getInputSignalName(int feature) override;

    private:
        int howMany;
        int wordSize;
	};


}//namespace

#endif
