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

#include "NeuralNetworks/Layers/Layer.hpp"

namespace flopoco {

	// new operator class declaration
    class PoolingLayer : public Layer {

    public:
        PoolingLayer(Target* target, NeuralNetwork* parent, int wordSize_, int horizontalSize_, int verticalSize_, int numberOfOutputFeatures_, bool calcAllParallel_, int paddingTop_, string paddingType_, int windowSize_, string activationFunction_="ReLU", int stride_=1, int paddingBot_=-1, int paddingLeft_=-1, int paddingRight_=-1, bool inputMemoryParallelAccess_=true, bool outputMemoryParallelAccess_=true);
        PoolingLayer(Target* target, NeuralNetwork* parent, LayerArguments* args, bool inputMemoryParallelAccess_=true, bool outputMemoryParallelAccess_=true);

		// destructor
        ~PoolingLayer() {}

        int getNumberOfInstances();

        virtual string getOutputSignalName(int feature) override;
        virtual string getInputSignalName(int feature) override;

    protected:
        virtual void generateVHDLCode(Target* target) override;

    private:
    	void createInputMemoryReset(Target* target);
	};


}//namespace

#endif
