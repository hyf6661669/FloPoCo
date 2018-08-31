#ifndef LAYERARGUMENTS_H
#define LAYERARGUMENTS_H

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

    class LayerArguments {

    public:
        LayerArguments();
        LayerArguments(string layerType_, int coreSize_, int inputHeight_, int inputWidth_, int inputDepth_, int wordSize_, int fraction_, int weightWordSize_, int weightFraction_, int numberOfOutputFeatures_, vector<double> weights_, int padding_, string paddingType_, bool inputFeaturesParallel_, bool outputFeaturesParallel_, string activationFunction_, int stride_, string id_);


		// destructor
        ~LayerArguments() {}

        void printLayerArguments();

        string getLayerType();
        int getCoreSize();
        int getInputHeight();
        int getInputWidth();
        int getInputDepth();
        int getWordSize();
        int getFraction();
        int getWeightWordSize();
        int getWeightFraction();
        int getNumberOfOutputFeatures();
        vector<double> getWeights();
        double getSpecificWeight(unsigned int index);
        double getConvWeight(unsigned int inputFeature, unsigned int outputFeature, unsigned int index);
        vector <vector <vector <double>>> getConvWeights();
        int getPaddingTop();
        int getPaddingBot();
        int getPaddingLeft();
        int getPaddingRight();
        int getPadding();
        string getPaddingType();
		bool getInputFeaturesParallel();
		bool getOutputFeaturesParallel();
		string getActivationFunction();
        int getStride();
        string getId();

        void setLayerType(string lt);
        void setCoreSize(int cs);
        void setInputHeight(int ih);
        void setInputWidth(int iw);
        void setInputDepth(int id);
        void setWordSize(int ws);
        void setFraction(int f);
        void setNumberOfOutputFeatures(int o);
        void setWeightWordSize(int wws);
        void setWeightFraction(int wf);
        void setWeights(vector<double> w);
        void setConvWeights(vector<vector<vector<double>>> w);
        void addWeight(double w);
        void setPadding(int p);
        void setPaddingTop(int p);
        void setPaddingBot(int p);
        void setPaddingLeft(int p);
        void setPaddingRight(int p);
        void setPaddingType(string p);
        void setInputFeaturesParallel(bool i);
        void setOutputFeaturesParallel(bool o);
        void setActivationFunction(string a);
        void setStride(int s);
        void setId(string i);

    protected:
        string layerType; // can be: Convolutional or FullConnected or Pooling
        int coreSize; // relevant for Convolution and Pooling
		int inputHeight;
		int inputWidth;
        int inputDepth; // = number of inputFeatures
		int wordSize;
		int fraction;
		int weightWordSize;
		int weightFraction;
        int numberOfOutputFeatures;
        vector<double> weights; // if this layer is a convLayer this vector can be seen to have following structure: weights[inputFeature][outputFeature][weightIndex] (use getter to get specific weights individually)
        int paddingTop; // relevant for Convolution and Pooling
        int paddingBot; // relevant for Convolution and Pooling
        int paddingLeft; // relevant for Convolution and Pooling
        int paddingRight; // relevant for Convolution and Pooling
        string paddingType; // relevant for Convolution and Pooling, can be "Zero" or "Value"
		bool inputFeaturesParallel;
		bool outputFeaturesParallel;
		string activationFunction;
        int stride;
        string id;

	private:
			
	};


}//namespace

#endif
