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
        LayerArguments(string layerType_, int coreSize_, int inputHeight_, int inputWidth_, int inputDepth_, int wordSize_, int fraction_, int weightWordSize_, int weightFraction_, int numberOfOutputFeatures_, vector<double> weights_, vector<double> biases_, int padding_, string paddingType_, bool calcAllParallel_, string activationFunction_, int stride_, unsigned int startAddress_, string id_);


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
        double getFullConnectedWeight(unsigned int inputNumber, unsigned int neuronNumber);
        vector <vector <double>> getFullConnectedWeights();
        unsigned int getWeightVectorSize();
        vector<double> getBiases();
        double getBias(unsigned int index);
        int getPaddingTop();
        int getPaddingBot();
        int getPaddingLeft();
        int getPaddingRight();
        int getPadding();
        string getPaddingType();
        bool getCalcAllParallel();
		string getActivationFunction();
        int getStride();
        unsigned int getStartAddress();
        string getId();
        bool getLutBasedAddressCalculation();

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
        void makeWeightsVectorBigger(unsigned int i);
        void setBiases(vector<double> b);
        void addBias(double b);
        void setPadding(int p);
        void setPaddingTop(int p);
        void setPaddingBot(int p);
        void setPaddingLeft(int p);
        void setPaddingRight(int p);
        void setPaddingType(string p);
        void setCalcAllParallel(bool p);
        void setActivationFunction(string a);
        void setStride(int s);
		void setStartAddress(unsigned int a);
        void setId(string i);
        void setLutBasedAddressCalculation(bool b);

		void reorderFullConnectedWeightsAndBiases();

		static vector<double> flatten4DTensor(vector<vector<vector<vector<double>>>> t);
		static vector<double> flatten2DTensor(vector<vector<double>> t);

    protected:
		void reorderFullConnectedWeights();
		void reorderFullConnectedBiases();

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
        vector<double> biases;
        int paddingTop; // relevant for Convolution and Pooling
        int paddingBot; // relevant for Convolution and Pooling
        int paddingLeft; // relevant for Convolution and Pooling
        int paddingRight; // relevant for Convolution and Pooling
        string paddingType; // relevant for Convolution and Pooling, can be "Zero" or "Value"
        bool calcAllParallel;
		string activationFunction;
        int stride;
        unsigned int startAddress; // relevant for Convolution and FullConnected
		bool fullConnectedWeightsAndBiasesOrdered; // relevant for FullConnected
        string id;
        bool lutBasedAddressCalculation; // relevant if calcAllParallel==false

	private:
			
	};


}//namespace

#endif
