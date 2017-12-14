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

namespace flopoco {

	// new operator class declaration
    class ConvolutionalLayer : public Operator {

    public:
        ConvolutionalLayer(Target* target, unsigned int wordSize_, unsigned int fraction_, unsigned int horizontalSize_, unsigned int verticalSize_, unsigned int windowSize_, unsigned int numberOfFeatures_, unsigned int numberOfInputFeatures_, vector <vector <vector <double> > > weights_, unsigned int weightWordSize_, unsigned int weightFraction_, int paddingLeft_, string paddingType_, unsigned int stride_=1, bool useAdderTree_=true, bool roundAfterConvCore_=true, string roundingType_="Truncation", string activationFunction_="ReLU", int paddingRight_=-1, int paddingTop_=-1, int paddingBot_=-1, string id_="0");

		// destructor
        ~ConvolutionalLayer() {}

        unsigned int getNumberOfInputFeatures(){return this->numberOfInputFeatures;}
        unsigned int getNumberOfOutputFeatures(){return this->numberOfFeatures;}

    private:
        void buildShiftAndPad(Target* target);
        void buildConvCores(Target* target);
        void buildAdderTrees(Target* target);
        void roundOutput(unsigned int wordSizeFrom, unsigned int fractionFrom, string signalNameFrom, string signalNameTo, string round="Truncation");
        void buildActivationFunction(Target* target);

        unsigned int wordSize;
        unsigned int fraction;
        unsigned int horizontalSize;
        unsigned int verticalSize;
        unsigned int windowSize;
        unsigned int numberOfFeatures;
        unsigned int numberOfInputFeatures;
        vector <vector <vector <double> > > weights; // format: weights[inputFeatureNo.][outputFeatureNo.][weightNo.]
        unsigned int weightWordSize;
        unsigned int weightFraction;
        int paddingLeft;
        int paddingRight;
        int paddingTop;
        int paddingBot;
        string paddingType;
        unsigned int stride;
		string id;
        bool useAdderTree;
        bool roundAfterConvCore;
        string roundingType;
        string activationFunction;

        vector <vector <ConvolutionalCoreSimple*>> ConvCores; // format: ConvCores[inputFeatureNo.][outputFeatureNo.]
        vector <IntAdderTree*> AdderTrees; // format: AdderTrees[outputFeatureNo.]
    };


}//namespace

#endif
