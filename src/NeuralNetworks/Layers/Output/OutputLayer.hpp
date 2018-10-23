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
#include "FindIndexOfMaximum.hpp"

namespace flopoco {

	// new operator class declaration
    class OutputLayer : public Layer {

    public:
        OutputLayer(Target* target, NeuralNetwork* parent, int numberOfFeatures_, int wordSize_, bool bramOnOutput_=false, int width_=-1, int height_=-1, bool threeBRAMS_=true, bool needIntermediateValue_=false, bool hasSpecificAddressOnResetPort=false, bool onlyOutputClassIndex=false);

		// destructor
        ~OutputLayer() {}

        virtual string getOutputSignalName(int feature) override;
        virtual string getInputSignalName(int feature) override;

        inline bool hasIntermediateValuePorts() const {return this->needIntermediateValue;}
		unsigned int getIndexWordSize() const;

    private:
        int numberOfFeatures;
        int wordSize;
        bool bramOnOutput;
        int width;
        int height;
        bool threeBRAMS;
        bool needIntermediateValue;

		FindIndexOfMaximum* myMaxPtr;
	};


}//namespace

#endif
