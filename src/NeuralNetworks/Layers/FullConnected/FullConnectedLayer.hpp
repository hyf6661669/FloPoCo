#ifndef FULLCONNECTEDLAYER_H
#define FULLCONNECTEDLAYER_H
/* Each Operator declared within the flopoco framework has 
   to inherit the class Operator and overload some functions listed below*/
#include "NeuralNetworks/Layers/Layer.hpp"
#include "NeuralNetworks/LayerArguments.hpp"

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
	class FullConnectedLayer : public Layer {

    public:
        FullConnectedLayer(Target* target, NeuralNetwork* parent, LayerArguments* lA, unsigned int weightsPerDMAAccess_, char roundingType_=0x00);

		// destructor
        ~FullConnectedLayer() {}

        string getInputSignalName(int number);
        string getOutputSignalName(int number);

    protected:
        virtual void generateVHDLCode(Target* target) override;

	private:
        void buildInputsAndOutputs(Target* target);
        void declareSignals(Target* target);
        void buildBiasLut(Target* target);
        void buildFullConnectedCore(Target* target);
        void buildController(Target* target);
        void buildWeightFetcher(Target* target);
        void buildWeightRegisters(Target* target);
        void roundOutput(Target* target);
        void buildActivationFunction(Target* target);

        void buildFullConnectedCoreAdder(Target* target);

        void resizeBias();

        unsigned int weightsPerDMAAccess;
        char roundingType;

        unsigned int intermediateResultWidth;
        unsigned int neuronCounterWidth;
        unsigned int weightCounterWidth;
        unsigned int weightsPerNeuron;
        unsigned int fullConnectedCorePipelineDepth;
        int fullConnectedCoreAdderPipelineDepth;
	};


}//namespace

#endif
