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
#include "IntAddSubCmp/IntAdderTree.hpp"
#include "NeuralNetworks/Layers/Layer.hpp"

namespace flopoco {

	// new operator class declaration
    class ConvolutionalLayer : public Layer {

    public:
        ConvolutionalLayer(Target* target, int wordSize_, int fraction_, int horizontalSize_, int verticalSize_, int windowSize_, int numberOfOutputFeatures_, int numberOfInputFeatures_, vector <vector <vector <double> > > weights_, int weightWordSize_, int weightFraction_, int paddingLeft_, string paddingType_, bool inputFeaturesParallel_, bool outputFeaturesParallel_, string id_, int stride_=1, bool useAdderTree_=true, bool roundAfterConvCore_=true, string roundingType_="Truncation", string activationFunction_="ReLU", int paddingRight_=-1, int paddingTop_=-1, int paddingBot_=-1);
        ConvolutionalLayer(Target* target, LayerArguments* la, bool useAdderTree_=true, bool roundAfterConvCore_=true, string roundingType_="Truncation");

		// destructor
        ~ConvolutionalLayer() {}

        //unsigned int getNumberOfInputFeatures(){return this->numberOfInputFeatures;}
        //unsigned int getNumberOfOutputFeatures(){return this->numberOfOutputFeatures;}

        virtual string getOutputSignalName(int feature) override;
        virtual string getInputSignalName(int feature) override;

    protected:
        virtual void generateVHDLCode(Target* target) override;
    private:
        void buildShiftAndPad(Target* target);
        void buildConvCores(Target* target);
        void buildAdderTrees(Target* target);
        void roundOutput(int wordSizeFrom, int fractionFrom, string signalNameFrom, string signalNameTo, string round="Truncation");
        void buildActivationFunction(Target* target);

        //unsigned int wordSize;
        //unsigned int fraction;
        //unsigned int horizontalSize;
        //unsigned int verticalSize;
        //unsigned int windowSize;
        //unsigned int numberOfOutputFeatures;
        //unsigned int numberOfInputFeatures;
        //vector <vector <vector <double> > > weights; // format: weights[inputFeatureNo.][outputFeatureNo.][weightNo.]
        //unsigned int weightWordSize;
        //unsigned int weightFraction;
        //int paddingLeft;
        //int paddingRight;
        //int paddingTop;
        //int paddingBot;
        //string paddingType;
        //unsigned int stride;
        //bool inputFeaturesParallel;
        //bool outputFeaturesParallel;
        //string id;
        bool useAdderTree;
        bool roundAfterConvCore;
        string roundingType;
        //string activationFunction;

        vector <vector <ConvolutionalCoreSimple*>> ConvCores; // format: ConvCores[inputFeatureNo.][outputFeatureNo.]
        vector <IntAdderTree*> AdderTrees; // format: AdderTrees[outputFeatureNo.]
    };


}//namespace

#endif
