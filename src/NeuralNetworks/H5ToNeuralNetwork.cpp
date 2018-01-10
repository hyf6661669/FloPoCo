#ifndef H5TONEURALNETWORK_H
#define H5TONEURALNETWORK_H
// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

// for reading the .txt file
#include <fstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// for random numbers
#include <cstdlib>
#include <ctime>

// include the header of the Operator
#include "H5ToNeuralNetwork.hpp"

// include the headers of all possible NN-Layers
#include "NeuralNetworks/Layers/FullConnected/FullConnectedLayer.hpp"
#include "NeuralNetworks/Layers/Convolutional/ConvolutionalLayer.hpp"
#include "NeuralNetworks/Layers/Convolutional/ConvolutionalCoreSimple.hpp"
#include "NeuralNetworks/ActivationFunctions/ReLU.hpp"
#include "NeuralNetworks/MemoryManagement/ControlledMemory.hpp"
#include "NeuralNetworks/Utility/PaddingGenerator.hpp"
#include "NeuralNetworks/Utility/GlobalController.hpp"
#include "NeuralNetworks/Layers/Pooling/PoolingLayer.hpp"

#include <ctime>


using namespace std;
namespace flopoco {




    H5ToNeuralNetwork::H5ToNeuralNetwork(Target* target, string pathToH5_) : Operator(target), pathToH5(pathToH5_) {

        // to calc needed time
        time_t timerStart = std::time(nullptr);

		// definition of the source file name, used for info and error reporting using REPORT 
		srcFileName="H5ToNeuralNetwork";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

		// definition of the name of the operator
		ostringstream name;
		name << "H5ToNeuralNetwork";
        setName(name.str());

        cout << "### pathToH5=" << pathToH5 << endl;
        ifstream txtfile;
        txtfile.open(pathToH5.c_str());

        if(!txtfile)
        {
            cout << "### Could not open the file. Your path may be wrong..." << endl;
        }
        else
        {
            cout << "### Yay! The File is open! :)" << endl;
        }

        // read Neural Network
        map <string, LayerArguments*> layers;
        vector <pair<string, string>> edges;

        // Input Layer
        LayerArguments* lAI = new LayerArguments();
        lAI->setLayerType("Input");
        layers["InputLayer"] = lAI;

        // ConvLayer 1
        LayerArguments* lA1 = new LayerArguments();
        lA1->setLayerType("Convolutional");
        lA1->setWordSize(8);
        lA1->setFraction(4);
        lA1->setInputDepth(1);
        lA1->setInputHeight(20);
        lA1->setInputWidth(40);
        lA1->setNumberOfOutputFeatures(30);
        lA1->setWeightWordSize(8);
        lA1->setWeightFraction(4);
        lA1->setInputFeaturesParallel(true);
        lA1->setOutputFeaturesParallel(true);
        lA1->setCoreSize(5);
        lA1->setConvWeights(H5ToNeuralNetwork::getRandomWeights(lA1->getInputDepth(),lA1->getNumberOfOutputFeatures(),lA1->getCoreSize()*lA1->getCoreSize(),0,(1<<lA1->getWeightWordSize()),lA1->getWeightFraction()));
        lA1->setPadding(2);
        lA1->setPaddingType("Value");
        lA1->setStride(1);
        lA1->setId("1");
        lA1->setActivationFunction("ReLU");
        layers["ConvLayer1"] = lA1;
        layers["ConvLayer1"]->printLayerArguments();

        // ConvLayer 2
        LayerArguments* lA2 = new LayerArguments();
        lA2->setLayerType("Convolutional");
        lA2->setWordSize(8);
        lA2->setFraction(4);
        lA2->setInputDepth(30);
        lA2->setInputHeight(20);
        lA2->setInputWidth(40);
        lA2->setNumberOfOutputFeatures(10);
        lA2->setWeightWordSize(8);
        lA2->setWeightFraction(4);
        lA2->setInputFeaturesParallel(true);
        lA2->setOutputFeaturesParallel(true);
        lA2->setCoreSize(3);
        lA2->setConvWeights(H5ToNeuralNetwork::getRandomWeights(lA2->getInputDepth(),lA2->getNumberOfOutputFeatures(),lA2->getCoreSize()*lA2->getCoreSize(),0,(1<<lA2->getWeightWordSize()),lA2->getWeightFraction()));
        lA2->setPadding(1);
        lA2->setPaddingType("Value");
        lA2->setStride(1);
        lA2->setId("2");
        lA2->setActivationFunction("ReLU");
        layers["ConvLayer2"] = lA2;
        layers["ConvLayer2"]->printLayerArguments();

        // Pooling Layer 1
        LayerArguments* lA3 = new LayerArguments();
        lA3->setLayerType("Pooling");
        lA3->setWordSize(8);
        lA3->setFraction(4);
        lA3->setInputDepth(10);
        lA3->setInputHeight(20);
        lA3->setInputWidth(40);
        lA3->setNumberOfOutputFeatures(10);
        //lA3->setWeightWordSize(8);
        //lA3->setWeightFraction(4);
        lA3->setInputFeaturesParallel(true);
        lA3->setOutputFeaturesParallel(true);
        lA3->setCoreSize(2);
        //lA3->setConvWeights(H5ToNeuralNetwork::getRandomWeights(lA->getInputDepth(),lA->getNumberOfOutputFeatures(),lA->getCoreSize()*lA->getCoreSize(),0,(1<<lA->getWeightWordSize()),lA->getWeightFraction()));
        lA3->setPaddingLeft(1);
        lA3->setPaddingTop(1);
        lA3->setPaddingType("Value");
        lA3->setStride(2);
        lA3->setId("1");
        layers["PoolingLayer1"] = lA3;
        layers["PoolingLayer1"]->printLayerArguments();

        // Output Layer
        LayerArguments* lAO = new LayerArguments();
        lAO->setLayerType("Output");
        layers["OutputLayer"] = lAO;

        // connections
        pair<string, string> edge;
        edge.first = "InputLayer";
        edge.second = "ConvLayer1";
        edges.push_back(edge);

        edge.first = edge.second;
        edge.second = "ConvLayer2";
        edges.push_back(edge);

        edge.first = edge.second;
        edge.second = "PoolingLayer1";
        edges.push_back(edge);

        edge.first = edge.second;
        edge.second = "OutputLayer";
        edges.push_back(edge);

        cout << "###    GOING INTO NEURAL NETWORK CONSTRUCTOR NOW!" << endl;
        // Constructor
        NeuralNetwork* NN = new NeuralNetwork(target,layers,edges,"TestNetwork");
        this->addSubComponent(NN);


        cout << "###    TIME NEEDED TO EXECUTE: " << ((float)std::difftime(std::time(nullptr),timerStart))/60 << " MINUTES!" << endl;

    }


    OperatorPtr H5ToNeuralNetwork::parseArguments(Target *target, vector<string> &args) {
        string pathToH5;
        UserInterface::parseString(args, "pathToH5", &pathToH5);
        return new H5ToNeuralNetwork(target, pathToH5);
    }
	
    void H5ToNeuralNetwork::registerFactory(){
        UserInterface::add(
                    "H5ToNeuralNetwork", // name
                    "Generate VHDL-Code from the .txt dump of a HDF5-File which describes a Neural Network (can be obtained via 'hdf5dump YourHDF5File.h5')", // description, string
                    "NeuralNetworks", // category, from the list defined in UserInterface.cpp
                    "", //seeAlso
                    // Now comes the parameter description string.
                    // Respect its syntax because it will be used to generate the parser and the docs
                    // Syntax is: a semicolon-separated list of parameterDescription;
                    // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
                    "pathToH5(string): the path to your .txt dump of the HDF5-file",
                    "",
                    H5ToNeuralNetwork::parseArguments
                    );
    }

    vector<vector<vector<double> > > H5ToNeuralNetwork::getRandomWeights(int a, int b, int c, int min, int max, int fractionalPart)
    {
        if(max-min<=0)
        {
            stringstream e;
            e << "Max<Min!";
            THROWERROR(e);
        }
        if(fractionalPart>=31)
        {
            stringstream e;
            e << "Fractional Part too big!";
            THROWERROR(e);
        }
        vector<vector<vector<double>>> level3;
        srand(time(0));
        for(int h=0;h<a;h++)
        {
            vector<vector<double>> level2;
            for(int w=0;w<b;w++)
            {
                vector<double> level1;
                for(int d=0;d<c;d++)
                {
                    int upperBound=rand()%max;
                    int lowerBound;
                    if(min>0)
                    {
                        lowerBound=rand()%min;
                    }
                    else
                    {
                        min=0;
                    }
                    int randNum=upperBound-lowerBound;
                    int upperBoundForScalingFactor=(1<<fractionalPart)-1;
                    int scalingFactor = rand()%upperBoundForScalingFactor+1;
                    float weight = ((float)randNum)/((float)scalingFactor);
                    level1.push_back(weight);
                }
                level2.push_back(level1);
            }
            level3.push_back(level2);
        }
        return level3;
    }

}//namespace flopoco
#endif
