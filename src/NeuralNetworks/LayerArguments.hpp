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
        LayerArguments(string layerType_, int number_, int coreSize_, int inputHeight_, int inputWidth_, int inputDepth_, int wordSize_, int fraction_, int weightWordSize_, int weightFraction_, int numberOfOutputFeatures_, vector<double> weights_, string id_);

		// destructor
        ~LayerArguments() {}

        string getLayerType();
        int getNumber();
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
        string getId();

        void setLayerType(string lt);
        void setNumber(int n);
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
        void setId(string i);
    private:
        string layerType; // can be: Convolutional or FullConnected or Pooling
		int number;
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
        string id;

	private:
			
	};


}//namespace

#endif
