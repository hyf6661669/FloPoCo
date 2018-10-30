#ifndef CONVOLUTIONALLAYER_H
#define CONVOLUTIONALLAYER_H
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

#include "ConvolutionalCoreSimple.hpp"
#include "NeuralNetworks/Layers/Convolutional/ConvolutionalCoreWithWeightInputs.hpp"
#include "IntAddSubCmp/IntAdderTree.hpp"
#include "NeuralNetworks/Layers/Layer.hpp"

namespace flopoco {

	// new operator class declaration
    class ConvolutionalLayer : public Layer {

    public:
        ConvolutionalLayer(Target* target, NeuralNetwork* parent, int wordSize_, int fraction_, int horizontalSize_, int verticalSize_, int windowSize_, int numberOfOutputFeatures_, int numberOfInputFeatures_, vector <vector <vector <double> > > weights_, vector<double> biases_, int weightWordSize_, int weightFraction_, int paddingLeft_, string paddingType_, bool calcAllParallel_, string id_, int stride_=1, bool useAdderTree_=true, bool roundAfterConvCore_=true, char roundingType_=0x00, string activationFunction_="ReLU", int paddingRight_=-1, int paddingTop_=-1, int paddingBot_=-1, bool inputMemoryParallelAccess_=true, bool outputMemoryParallelAccess_=true);
        ConvolutionalLayer(Target* target, NeuralNetwork* parent, LayerArguments* la, bool useAdderTree_=true, bool roundAfterConvCore_=true, char roundingType_=0x00, bool inputMemoryParallelAccess_=true, bool outputMemoryParallelAccess_=true);

		// destructor
        ~ConvolutionalLayer() {}

        virtual string getOutputSignalName(int feature) override;
        virtual string getInputSignalName(int feature) override;
        virtual string getIntermediateResultName(unsigned int featureNumber) override;

        int getInputCounterWidth() const;
        int getOutputCounterWidth() const;

        static string getInputFeatureCounterOutputSignalName() {return "InputFeatureCounter_out";}
        static string getOutputFeatureCounterOutputSignalName() {return "OutputFeatureCounter_out";}
        static string getConvCoreInputSignalName(unsigned int inputPort) {return "ConvCoreDataIn_"+to_string(inputPort);}
        static string getConvCoreWeightSignalName(unsigned int inputPort) {return "ConvCoreWeight_"+to_string(inputPort);}
        static string getConvCoreOutputSignalName() {return "ConvCoreDataOut";}

    protected:
        virtual void generateVHDLCode(Target* target) override;
    private:
        void buildShiftAndPadParallel(Target* target);
        void buildConvCoresParallel(Target* target);
        void buildAdderTreesParallel(Target* target);
        void roundOutput(Target* target, int wordSizeFrom, int fractionFrom, string signalNameFrom, string signalNameTo, string instanceName, char round=0x00);
        void buildActivationFunction(Target* target, int outFeatureCounter=0);

        void buildParallelCalculation(Target* target);

        void buildSerialCalculation(Target* target);
        void createPortsSerial(Target* target);
        void declareSignals(Target* target);
        void createDataGuardsAndMultiplexersSerial(Target* target);
        ConvolutionalCoreWithWeightInputs* createConvCoreSerial(Target* target);
        void createBiasesSerial(Target* target);
        unsigned int createBitHeapSerial(Target* target, unsigned int wIn);
        void createPipelineRegistersSerial(Target* target, unsigned int pipelineStagesConv, unsigned int pipelineStagesBitH);
        void createShitAndPadSerial(Target* target);
        void createDMAControllerSerial(Target* target);
        void createCountersSerial(Target* target);
        void createActivationFunctionSerial(Target* target);
        void createLutForIntermediateReset(Target* target);

        void createGetNewDataAddressLut(Target* target);

        bool useAdderTree;
        bool roundAfterConvCore;
        char roundingType;

        vector <vector <ConvolutionalCoreSimple*>> ConvCores; // format: ConvCores[inputFeatureNo.][outputFeatureNo.]
        vector <IntAdderTree*> AdderTrees; // format: AdderTrees[outputFeatureNo.]
    };


}//namespace

#endif
