#ifndef H5TONEURALNETWORK_H
#define H5TONEURALNETWORK_H
// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

//for reading the .txt file
#include <fstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "H5ToNeuralNetwork.hpp"

// include the headers of all possible NN-Layers
#include "FullConnectedLayer.hpp"
#include "ConvolutionalLayer.hpp"
#include "ConvolutionalCoreSimple.hpp"
#include "ReLU.hpp"
#include "ControlledMemory.hpp"
#include "PaddingGenerator.hpp"
#include "GlobalController.hpp"
#include "PoolingLayer.hpp"

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

        cout << "pathToH5=" << pathToH5 << endl;
        ifstream txtfile;
        txtfile.open(pathToH5.c_str());

        if(!txtfile)
        {
            cout << "Could not open the file. Your path may be wrong..." << endl;
        }
        else
        {
            cout << "Yay! The File is open! :)" << endl;
        }

        //parameters for a ConvolutionalCoreSimple
        unsigned int wordSize=8;
        int fraction=8;
        unsigned int weightWordSize=8;
        int weightFraction=8;
        unsigned int size=5;
        vector<double> weights;
        weights.push_back(1.25);
        weights.push_back(2.0);
        weights.push_back(3.5);
        weights.push_back(4.0);
        weights.push_back(5.0);
        weights.push_back(6.125);
        weights.push_back(7.625);
        weights.push_back(8.875);
        weights.push_back(9.0);
        weights.push_back(1.25);
        weights.push_back(2.0);
        weights.push_back(3.5);
        weights.push_back(4.0);
        weights.push_back(5.0);
        weights.push_back(6.125);
        weights.push_back(7.625);
        weights.push_back(8.875);
        weights.push_back(9.0);
        weights.push_back(1.25);
        weights.push_back(2.0);
        weights.push_back(3.5);
        weights.push_back(4.0);
        weights.push_back(5.0);
        weights.push_back(6.125);
        weights.push_back(7.625);

        //ports and signals in Toplevel
        for(unsigned int i=0; i<size*size; i++)
        {
            this->addInput("i"+to_string(i), wordSize);
            this->vhdl << tab << declare("i"+to_string(i)+"_out", wordSize) << " <= i" << i << ";" << endl;
        }
        this->addOutput("o", wordSize);

        //building testOp
        cout << "ConvCoreSimple \n";
        Operator* testOp = new ConvolutionalCoreSimple(target, wordSize, fraction, weightWordSize, weightFraction, size, weights,true,"Truncation","testoMesto");

        //adding component
        this->addSubComponent(testOp);

        //inPortMap
        for(unsigned int i=0; i<size*size; i++)
        {
            string portInOp, actualSignalName;
            portInOp = "X" + to_string(i);
            actualSignalName = "i" + to_string(i) + "_out";
            inPortMap(testOp, portInOp, actualSignalName);
        }
//        addInput("ConvCoreSimple_validData_i",1);
//        inPortMap(testOp,"validData_i","ConvCoreSimple_validData_i");
//        addOutput("ConvCoreSimple_validData_o",1);
//        outPortMap(testOp,"validData_o","ConvCoreSimple_validData_o",false);

        //outPortMap
        declare("ConvCore_Out", wordSize);
        outPortMap(testOp, "R", "ConvCore_out", false);

        //instance
        this->vhdl << this->instance(testOp, "ConvCore_instance");

        //ReLU
        Operator* testOp2 = new ReLU(target, wordSize);
        this->addSubComponent(testOp2);
        inPortMap(testOp2, "X", "ConvCore_Out");
        declare("ReLU_out", wordSize);
        outPortMap(testOp2, "R", "ReLU_out", false);
        vhdl << instance(testOp2, "ReLU_instance");

        //outport <= ReLU_out
        vhdl << "o <= ReLU_out;" << endl;

        //ControlledMemory
        unsigned int wAddr=13;
        unsigned int dataWidth=16;
        cout << "ControlledMemory \n";
        Operator* testOp3 = new ControlledMemory(target, dataWidth, wAddr);
        this->addSubComponent(testOp3);

        addInput("ControlledMemoryNewStep",1);
        inPortMap(testOp3, "newStep", "ControlledMemoryNewStep");
        addInput("ControlledMemoryValidData_i",1);
        inPortMap(testOp3, "validData_i", "ControlledMemoryValidData_i");
        addInput("ControlledMemoryGetNewData",1);
        inPortMap(testOp3, "getNewData", "ControlledMemoryGetNewData");
        addInput("ControlledMemoryData_i",dataWidth);
        inPortMap(testOp3, "data_i", "ControlledMemoryData_i");

        addOutput("ControlledMemoryValidData_o",1);
        outPortMap(testOp3, "validData_o", "ControlledMemoryValidData_o", false);
        addOutput("ControlledMemoryData_o",dataWidth);
        outPortMap(testOp3, "data_o", "ControlledMemoryData_o", false);

        vhdl << instance(testOp3, "ControlledMemory_instance");

        unsigned int padWordSize=16;
        unsigned int padWindowSize=7;
        unsigned int padHSize=30;
        unsigned int padVSize=23;
        unsigned int padStride=2;
        int padTop=2;
        int padBot=2;
        int padLeft=2;
        int padRight=2;
        string padType="Value";
        cout << "PaddingGenerator \n";
        PaddingGenerator* padOp = new PaddingGenerator(target, padWordSize, padWindowSize, padHSize, padVSize, padTop, padStride, padType, padBot, padLeft, padRight, true);
        this->addSubComponent(padOp);

        for(unsigned int i=0; i<padWindowSize; i++)
        {
            for(unsigned int j=0; j<padWindowSize; j++)
            {
                cout << "padOp " << "i" << i << "j" << j << endl;
                this->addInput("padX"+to_string(i)+to_string(j),padWordSize);
                this->inPortMap(padOp, padOp->inputNames[i][j], "padX"+to_string(i)+to_string(j));
                this->addOutput("padY"+to_string(i)+to_string(j),padWordSize);
                this->outPortMap(padOp, padOp->outputNames[i][j], "padY"+to_string(i)+to_string(j), false);
            }
        }
        this->addInput("padNewStep",1);
        this->addInput("padValidData_i",1);
        this->addOutput("padValidData_o",1);
        this->addOutput("padGetNewData",1);
        this->addOutput("padFinished",1);

        this->inPortMap(padOp,"newStep","padNewStep");
        this->inPortMap(padOp,"validData_i","padValidData_i");
        this->outPortMap(padOp,"validData_o","padValidData_o",false);
        this->outPortMap(padOp,"getNewData","padGetNewData",false);
        this->outPortMap(padOp,"finished","padFinished",false);

        this->vhdl << instance(padOp, "PaddingGenerator_instance");

        // global controller
        cout << "GlobalController \n";
        GlobalController* globalOp = new GlobalController(target, 5);
        this->addSubComponent(globalOp);
        for(unsigned int i=0; i<5; i++)
        {
            this->addInput("globalOpFinished_"+to_string(i),1);
            this->inPortMap(globalOp, "finished_"+to_string(i), "globalOpFinished_"+to_string(i));
        }
        this->addOutput("globalOpNewStep",1);
        this->outPortMap(globalOp,"newStep","globalOpNewStep",false);
        this->vhdl << instance(globalOp,"globalController_instance");


        // Convolutional Layer
        /**
        cout << "ConvLayer \n";
        unsigned int convLayerWordSize=16;
        unsigned int convLayerFraction=8;
        unsigned int convLayerHSize=30;
        unsigned int convLayerVSize=20;
        unsigned int convLayerWindowSize=3;
        unsigned int convLayerNumberOfFeatures=10;
        unsigned int convLayerNumberOfInputFeatures=5;
        vector<vector<vector<double>>> convLayerWeights;
        for(unsigned int inputFCounter=0;inputFCounter<convLayerNumberOfInputFeatures;inputFCounter++)
        {
            vector<vector<double>> tmp1;
            convLayerWeights.push_back(tmp1);
            for(unsigned int outputFCounter=0; outputFCounter<convLayerNumberOfFeatures; outputFCounter++)
            {
                vector<double>tmp2;
                for(unsigned int indexCounter=0; indexCounter<convLayerWindowSize*convLayerWindowSize;indexCounter++)
                {
                    //fill tmp2
                    tmp2.push_back(((double)inputFCounter-2)/((double)outputFCounter+1));
                }
                convLayerWeights[inputFCounter].push_back(tmp2);
            }
        }
        cout << "   weights (before constructor): " << endl;
        for(unsigned int inputFCounter=0;inputFCounter<convLayerWeights.size();inputFCounter++)
        {
            cout << "inputFeature: " << inputFCounter << endl;
            for(unsigned int outputFCounter=0;outputFCounter<convLayerWeights[inputFCounter].size();outputFCounter++)
            {
                cout << "   outputFeature: " << outputFCounter << endl;
                for(unsigned int indexCounter=0; indexCounter<convLayerWeights[inputFCounter][outputFCounter].size(); indexCounter++)
                {
                    cout << "       " << convLayerWeights[inputFCounter][outputFCounter][indexCounter] << endl;
                }
            }
            cout << endl;
        }
        unsigned int convLayerWeightWordSize=16;
        unsigned int convLayerWeightFraction=8;
        int convLayerPaddingLeft=1;
        int convLayerPaddingRight=1;
        int convLayerPaddingTop=1;
        int convLayerPaddingBot=1;
        string convLayerPaddingType="Value";
        ConvolutionalLayer* convLayer = new ConvolutionalLayer(target,convLayerWordSize,convLayerFraction,convLayerHSize,convLayerVSize,convLayerWindowSize,
                                                               convLayerNumberOfFeatures,convLayerNumberOfInputFeatures,convLayerWeights,convLayerWeightWordSize,convLayerWeightFraction,
                                                               convLayerPaddingLeft,convLayerPaddingType,2,true,true,"Truncation","ReLU",
                                                               convLayerPaddingRight,convLayerPaddingTop,convLayerPaddingBot,"0");
        this->addSubComponent(convLayer);

        this->addInput("convLayerNewStep",1);
        this->inPortMap(convLayer,"newStep","convLayerNewStep");
        for(unsigned int i=0; i<convLayerNumberOfInputFeatures;i++)
        {
            this->addInput("convLayerValidData_i_"+to_string(i),1);
            this->inPortMap(convLayer,"validData_i_"+to_string(i),"convLayerValidData_i_"+to_string(i));
            this->addInput("convLayerX"+to_string(i),convLayerWordSize);
            this->inPortMap(convLayer,"X"+to_string(i),"convLayerX"+to_string(i));
            this->addOutput("convLayerGetNewData"+to_string(i),1);
            this->outPortMap(convLayer,"getNewData"+to_string(i),"convLayerGetNewData"+to_string(i),false);
        }
        this->addOutput("convLayerFinished",1);
        this->outPortMap(convLayer,"finished","convLayerFinished",false);
        for(unsigned int i=0; i<convLayerNumberOfFeatures; i++)
        {
            this->addOutput("convLayerValidData_o_"+to_string(i),1);
            this->outPortMap(convLayer,"validData_o_"+to_string(i),"convLayerValidData_o_"+to_string(i),false);
            this->addOutput("convLayerR"+to_string(i),convLayerWordSize);
            this->outPortMap(convLayer,"R"+to_string(i),"convLayerR"+to_string(i),false);
        }
        this->vhdl << instance(convLayer, "CONVLAYER_SHITYOURPANTS");
*/
        // PoolingLayer
        unsigned int pooWordSize=16;
        unsigned int pooHSize=15;
        unsigned int pooVSize=11;
        unsigned int pooNumberOfFeatures=7;
        bool pooParallel=true;
        int pooPaddingTop=1;
        unsigned int pooWindowSize=2;
        unsigned int pooStride=2;
        int pooPaddingBot=0;
        int pooPaddingLeft=1;
        int pooPaddingRight=0;
        string pooPaddingType="Value";
cout << "### POOLING!" << endl;
        PoolingLayer* pooOp = new PoolingLayer(target,pooWordSize,pooHSize,pooVSize,pooNumberOfFeatures,pooParallel,pooPaddingTop,pooWindowSize,pooStride,pooPaddingBot,pooPaddingLeft,pooPaddingRight,pooPaddingType);
        addSubComponent(pooOp);
cout << "###    1" << endl;
        this->addInput("pooNewStep",1);
        inPortMap(pooOp,"newStep","pooNewStep");
        this->addOutput("pooFinished",1);
        outPortMap(pooOp,"finished","pooFinished",false);
cout << "###    2" << endl;
        for(unsigned int featC=0; featC<pooOp->getNumberOfInstances(); featC++)
        {
            cout << "###    3.1: " << featC << endl;
            addInput("pooX"+to_string(featC),pooWordSize);
            inPortMap(pooOp,"X_"+to_string(featC),"pooX"+to_string(featC));
            addOutput("pooR"+to_string(featC),pooWordSize);
            outPortMap(pooOp,"R_"+to_string(featC),"pooR"+to_string(featC),false);

            addInput("pooValidData_i"+to_string(featC),1);
            inPortMap(pooOp,"validData_i"+to_string(featC),"pooValidData_i"+to_string(featC));

            addOutput("pooGetNewData"+to_string(featC),1);
            outPortMap(pooOp,"getNewData"+to_string(featC),"pooGetNewData"+to_string(featC),false);
            addOutput("pooValidData_o"+to_string(featC),1);
            outPortMap(pooOp,"validData_o"+to_string(featC),"pooValidData_o"+to_string(featC),false);
            cout << "###    3.2: " << featC << endl;
        }
cout << "###    4" << endl;
        vhdl << instance(pooOp,"poolingShit") << endl;

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

}//namespace flopoco
#endif
