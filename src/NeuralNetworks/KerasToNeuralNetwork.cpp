#ifndef KERASTONEURALNETWORK_H
#define KERASTONEURALNETWORK_H
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
#include "KerasToNeuralNetwork.hpp"

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




    KerasToNeuralNetwork::KerasToNeuralNetwork(Target* target, string pathToFolder_, bool serial, int wordSize_, int fraction_, int weightWordSize_, int weightFraction_, char roundingType, bool onlyOutputClassIndex_) :
        Operator(target), pathToFolder(pathToFolder_), wordSize(wordSize_), fraction(fraction_), weightWordSize(weightWordSize_), weightFraction(weightFraction_), onlyOutputClassIndex(onlyOutputClassIndex_), orderedLayerNames() {

        // word sizes
        if(this->weightWordSize<0)
        {
            this->weightWordSize = this->wordSize;
        }
        if(this->weightFraction<0)
        {
            this->weightFraction = this->fraction;
        }

		// definition of the source file name, used for info and error reporting using REPORT 
		srcFileName="KerasToNeuralNetwork";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

		// definition of the name of the operator
		ostringstream name;
		name << "KerasToNeuralNetwork";
        setName(name.str());

        if(roundingType>1 || roundingType<0) THROWERROR("Undefined Rounding Type, only 0 (truncation) and 1 (saturation) are possible yet");

        NeuralNetwork* NN = this->buildFromCSVFiles(target,serial,roundingType);

        this->addSubComponent(NN);

    }


    OperatorPtr KerasToNeuralNetwork::parseArguments(Target *target, vector<string> &args) {
        string pathToFolder_;
        int wordSize_;
        int fraction_;
        int weigtWordSize_;
        int weightFraction_;
        bool serial_;
        int roundingType;
        bool onlyOutputClassIndex;
        UserInterface::parseString(args, "pathToFolder", &pathToFolder_, false);
        UserInterface::parsePositiveInt(args, "wordSize", &wordSize_, false);
        UserInterface::parsePositiveInt(args, "fraction", &fraction_, false);
        UserInterface::parseInt(args, "weightWordSize", &weigtWordSize_, false);
        UserInterface::parseInt(args, "weightFraction", &weightFraction_, false);
        UserInterface::parseBoolean(args, "serial", &serial_, false);
        UserInterface::parseInt(args, "roundingType", &roundingType, false);
        UserInterface::parseBoolean(args, "onlyOutputClassIndex", &onlyOutputClassIndex, false);
        return new KerasToNeuralNetwork(target, pathToFolder_, serial_, wordSize_, fraction_, weigtWordSize_, weightFraction_, (char)roundingType, onlyOutputClassIndex);
    }
	
    void KerasToNeuralNetwork::registerFactory(){
        UserInterface::add(
                    "KerasToNeuralNetwork", // name
                    "Generate VHDL-Code from .csv files and a toplevel .txt file, that describe your Neural Network", // description, string
                    "NeuralNetworks", // category, from the list defined in UserInterface.cpp
                    "", //seeAlso
                    // Now comes the parameter description string.
                    // Respect its syntax because it will be used to generate the parser and the docs
                    // Syntax is: a semicolon-separated list of parameterDescription;
                    // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
                    "pathToFolder(string): the path to your folder that contains files which describe your Neural Network; \
                    serial(bool)=true: the neural network uses time multiplexing to save MANY resources (must be TRUE if it contains a full connected layer!);\
                    wordSize(int)=16: input word size;\
                    fraction(int)=0: number of fraction-bits;\
                    weightWordSize(int)=-1: word size of the weights, negative number: same as word size;\
                    weightFraction(int)=-1: word size of the weights, negative number: same as fraction;\
                    roundingType(int)=0: 0 => truncation, 1 => saturation;\
                    onlyOutputClassIndex(bool)=false: true => do NOT create output-registers for the output neurons (only relevant if the last layer is a full connected layer)",
                    "",
                    KerasToNeuralNetwork::parseArguments
                    );
    }

    vector<vector<vector<double> > > KerasToNeuralNetwork::getRandomWeights(int a, int b, int c, int wordSize, int fractionalPart)
    {
        if(wordSize>=31 || wordSize<0)
        {
            stringstream e;
            e << "Invalid Word Size";
            THROWERROR(e.str());
        }
        if(fractionalPart>=31 || fractionalPart<0)
        {
            stringstream e;
            e << "Invalid Fractional Part!";
            THROWERROR(e.str());
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
                    int range=(1<<wordSize);
                    int randNum=rand()%range; // interval: [0, 2^wordSize)
                    int splitRange=(1<<(wordSize-1));
                    randNum=randNum-splitRange; // interval: [2^(wordSize-1), 2^(wordSize-1))
                    float weight = ((float)randNum)/((float)(1<<fractionalPart)); // interval: [2^(wordSize-1)*2^fractionalPart, 2^(wordSize-1)*2^fractionalPart)
                    level1.push_back(weight);
                    float min = -1*((float)(1<<(wordSize-1))/((float)(1<<fractionalPart)));
                    float max = ((float)((1<<(wordSize-1))-1))/((float)(1<<fractionalPart));
                }
                level2.push_back(level1);
            }
            level3.push_back(level2);
        }
        return level3;
    }

    NeuralNetwork* KerasToNeuralNetwork::buildTestForZedboard(Target* target)
    {
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
        lA1->setWordSize(9);
        lA1->setFraction(0);
        lA1->setInputDepth(1);
        lA1->setInputHeight(100);
        lA1->setInputWidth(100);
        lA1->setNumberOfOutputFeatures(1);
        lA1->setWeightWordSize(9);
        lA1->setWeightFraction(0);
        lA1->setCalcAllParallel(true);
        lA1->setCoreSize(3);
        // WEIGHTS START
        //lA1->setConvWeights(KerasToNeuralNetwork::getRandomWeights(lA1->getInputDepth(),lA1->getNumberOfOutputFeatures(),lA1->getCoreSize()*lA1->getCoreSize(),lA1->getWeightWordSize(),lA1->getWeightFraction()));
        vector<vector<vector<double>>> v;
        vector<vector<double>> v1;
        vector<double> v2;
        v2.push_back(-1.0);
        v2.push_back(-1.0);
        v2.push_back(-1.0);
        v2.push_back(-1.0);
        v2.push_back(8.0);
        v2.push_back(-1.0);
        v2.push_back(-1.0);
        v2.push_back(-1.0);
        v2.push_back(-1.0);
        v1.push_back(v2);
        v.push_back(v1);
        lA1->setConvWeights((v));
        // WEIGHTS END
        lA1->setPadding(1);
        lA1->setPaddingType("Value");
        lA1->setStride(1);
        lA1->setId("1");
        lA1->setActivationFunction("relu");
        layers["ConvLayer1"] = lA1;
        layers["ConvLayer1"]->printLayerArguments();

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
        edge.second = "OutputLayer";
        edges.push_back(edge);

        // Constructor
        NeuralNetwork* NN = new NeuralNetwork(target,layers,edges,"TestNetwork");

        return NN;

    }

    NeuralNetwork *KerasToNeuralNetwork::buildOneLayerForDebugging(Target *target)
    {

        NeuralNetwork* NN = new NeuralNetwork(target,"DebugNetwork");

        map <string, LayerArguments*> layers;
        vector <pair<string, string>> edges;

        // Input Layer
        LayerArguments* lAI = new LayerArguments();
        lAI->setLayerType("Input");
        layers["InputLayer"] = lAI;

        // ConvLayer 1
        LayerArguments* lA1 = new LayerArguments();
        lA1->setLayerType("Convolutional");
        lA1->setWordSize(9);
        lA1->setFraction(0);
        lA1->setInputDepth(1);
        lA1->setInputHeight(100);
        lA1->setInputWidth(100);
        lA1->setNumberOfOutputFeatures(1);
        lA1->setWeightWordSize(9);
        lA1->setWeightFraction(0);
        lA1->setCalcAllParallel(true);
        lA1->setCoreSize(3);
        // WEIGHTS START
        vector<vector<vector<double>>> v;
        vector<vector<double>> v1;
        vector<double> v2;
        v2.push_back(-1.0);
        v2.push_back(-1.0);
        v2.push_back(-1.0);
        v2.push_back(-1.0);
        v2.push_back(8.0);
        v2.push_back(-1.0);
        v2.push_back(-1.0);
        v2.push_back(-1.0);
        v2.push_back(-1.0);
        v1.push_back(v2);
        v.push_back(v1);
        lA1->setConvWeights((v));
        // WEIGHTS END
        lA1->setPadding(1);
        lA1->setPaddingType("Value");
        lA1->setStride(1);
        lA1->setId("1");
        lA1->setActivationFunction("relu");
        layers["ConvLayer1"] = lA1;
        layers["ConvLayer1"]->printLayerArguments();

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
        edge.second = "OutputLayer";
        edges.push_back(edge);

        // setting layers/edges
        NN->setLayers(layers);
        NN->setEdges(edges);

        // build vhdl!
        NN->buildOnlyOneLayer(target,"ConvLayer1");
        return NN;
    }

    NeuralNetwork *KerasToNeuralNetwork::buildFromCSVFiles(Target *target, bool serial, char roundingType)
    {
        string pathToToplevelFile = pathToFolder + "NeuralNetwork.txt";
        this->processToplevelFile(pathToToplevelFile);
        this->layerIDCounter = 0;
        for(unsigned int i=0; i<this->orderedLayerNames.size(); i++)
        {
            if(i>0)
            {
                this->setUndefinedInputShapes(i);
            }

            string layerName = this->orderedLayerNames[i];
            if(layerName!="Input" && layerName!="Output")
            {
                this->readCSV(layerName);
            }
            this->layerIDCounter++;
        }
        for(auto it : this->layerArgs)
        {
            LayerArguments* lA = it.second;
            if(lA->getLayerType()=="FullConnected")
            {
                lA->reorderFullConnectedWeightsAndBiases();
            }
        }

        NeuralNetwork* returnMe = new NeuralNetwork(target,this->layerArgs,this->edges,this->networkName,serial,roundingType,0x00,this->onlyOutputClassIndex);
        return returnMe;
    }

    void KerasToNeuralNetwork::processToplevelFile(string pathToToplevelFile)
    {
        ifstream toplevelfile;
        toplevelfile.open(pathToToplevelFile.c_str());

        if(!toplevelfile)
        {
            cout << "### Could not open the file. Your path may be wrong..." << endl;
            stringstream e;
            e << "Textfile not found, path to folder: '" << this->pathToFolder << "'; you may have forgotten the last slash/backslash?! (Path MUST end with a '/')";
            THROWERROR(e.str());
        }
        else
        {
            cout << "### Yay! The File is open! :)" << endl;
        }

        while(!toplevelfile.eof())
        {
            string linebuffer;
            getline(toplevelfile,linebuffer);
            istringstream lineAsStream(linebuffer);
            // get first keyword; can be 'name', 'layer' or 'connection'
            if(!lineAsStream)
            {
                continue;
            }
            string word;
            lineAsStream>>word;
            if(word=="name")
            {
                lineAsStream>>this->networkName;
            }
            else if(word=="layer")
            {
                string name;
                string type;
                lineAsStream>>type;
                lineAsStream>>name;
                this->layerArgs[name] = new LayerArguments();
                this->layerArgs[name]->setLayerType(type);

                // the word sizes aren't given in keras
                this->layerArgs[name]->setWordSize(this->wordSize);
                this->layerArgs[name]->setWeightWordSize(this->weightWordSize);
                this->layerArgs[name]->setFraction(this->fraction);
                this->layerArgs[name]->setWeightFraction(this->weightFraction);

                //this->layerTypes[name]=type; // setLayerType
            }
            else if(word=="connection")
            {
                pair<string,string> connection;
                lineAsStream>>connection.first;
                lineAsStream>>connection.second;
                this->edges.push_back(connection);
            }
            else if(word=="")
            {
                continue;
            }
            else
            {
                stringstream e;
                e << "Unknown keyword '" << word << "'";
                THROWERROR(e.str());
            }
        }
        this->orderedLayerNames = this->orderLayerNames(this->edges);
    }

    void KerasToNeuralNetwork::readCSV(string layerName)
    {
        // file with information, weights, etc.
        string pathToFile = this->pathToFolder+layerName+".csv";
        //ifstream file;
        //file.open(pathToFile.c_str());
        //if(!file)
        //{
            //stringstream e;
            //e << "Could not find the .csv-file for layer '" << layerName << "'";
            //THROWERROR(e.str());
        //}

        // header with additional information about the layer
        this->readHeaderOfCSV(pathToFile,layerName);
        // set the layer ID
        this->layerArgs[layerName]->setId(to_string(this->layerIDCounter));
        // fill weights!
        this->readWeightsOfCSV(pathToFile,layerName);

        pathToFile = this->pathToFolder+layerName+"_biases.csv";
        this->readBiasOfCSV(pathToFile,layerName);

        //file.close();
    }

    void KerasToNeuralNetwork::readHeaderOfCSV(string pathToFile, string layerName)
    {
        ifstream file;
        file.open(pathToFile.c_str());
        if(!file)
        {
            stringstream e;
            e << "Could not find the .csv-file for layer '" << layerName << "'";
            THROWERROR(e.str());
        }

        // read header
        string line;
        getline(file,line);
        while(line[0]=='#' || line=="nan")
        {
            if(line=="nan")
            {
                // file contains exactly "nan" after header (happens for pooling layers)
                this->layerArgs[layerName]->setNumberOfOutputFeatures(this->layerArgs[layerName]->getInputDepth());
                return;
            }
            istringstream lineAsStream(line);
            string keyword;
            lineAsStream >> keyword; // now keyword="#", ignore that
            lineAsStream >> keyword;

            // fill this->layerArgs, and expect the pointer to have a valid value!
            if(keyword == "kernel_size")
            {
                string value;
                lineAsStream >> value;
                this->layerArgs[layerName]->setCoreSize(stod(value));
            }
            else if(keyword=="filters")
            {
                string value;
                lineAsStream >> value;
                this->layerArgs[layerName]->setNumberOfOutputFeatures(stoi(value));
            }
            else if(keyword=="strides")
            {
                string valueH;
                string valueV;
                // keras supports two stride-values (horizontal and vertical)
                lineAsStream >> valueH;
                lineAsStream >> valueV;
                this->layerArgs[layerName]->setStride(stoi(valueH));
            }
            else if(keyword=="activation")
            {
                string value;
                lineAsStream >> value;
                this->layerArgs[layerName]->setActivationFunction(value);
            }
            else if(keyword=="input_shape")
            {
                string inputDepth;
                string inputWidth;
                string inputHeight;
                lineAsStream >> inputDepth;
                lineAsStream >> inputWidth;
                lineAsStream >> inputHeight;
                if(inputDepth!="not_specified")
                {
                    this->layerArgs[layerName]->setInputDepth(stoi(inputDepth));
                }
                if(inputWidth!="not_specified")
                {
                    this->layerArgs[layerName]->setInputWidth(stoi(inputWidth));
                }
                if(inputHeight!="not_specified")
                {
                    this->layerArgs[layerName]->setInputHeight(stoi(inputHeight));
                }
                // in case it's not specified, leave it at default value (-1), so the neural network class can handle it
            }
            else if(keyword=="padding")
            {
                string value;
                lineAsStream >> value;
                if(value=="same")
                {
                    // calculate needed padding to maintain image size
                    int coreS = this->layerArgs[layerName]->getCoreSize();
                    if(coreS<1)
                    {
                        stringstream e;
                        e << "Need a valid core size to calculate needed padding; current core size: '" << coreS << "'";
                        THROWERROR(e.str());
                    }
                    if(coreS%2==1)
                    {
                        this->layerArgs[layerName]->setPadding((coreS-1)/2);
                    }
                    else
                    {
                        this->layerArgs[layerName]->setPaddingTop(coreS/2);
                        this->layerArgs[layerName]->setPaddingLeft(coreS/2);
                        this->layerArgs[layerName]->setPaddingBot((coreS)/2-1);
                        this->layerArgs[layerName]->setPaddingRight((coreS/2)-1);
                        // check, what happens when stride, etc...
                    }

                    // the is currently only zero padding in keras
                    this->layerArgs[layerName]->setPaddingType("Zero");
                }
                else
                {
                    this->layerArgs[layerName]->setPadding(0);
                }
            }
            else
            {
                stringstream e;
                e << "Unknown keyword in csv-header: '" << keyword << "'";
                THROWERROR(e.str());
            }
            getline(file,line); // read next line
        }

        file.close();
    }

    void KerasToNeuralNetwork::readWeightsOfCSV(string pathToFile, string layerName)
    {
        ifstream file;
        file.open(pathToFile.c_str());

        string linebuffer;
        getline(file,linebuffer);

        // we don't care about the header in this function...
        while(linebuffer[0]=='#')
        {
            getline(file,linebuffer);
        }

        // this shit is a pooling layer!
        if(linebuffer=="nan") return;

        LayerArguments* lA = this->layerArgs[layerName];
        int inputDepth = lA->getInputDepth();
        //int inputHeight = this->layerArgs[layerName]->getInputHeight();
        //int inputWidth = this->layerArgs[layerName]->getInputWidth();
        int outputDepth = lA->getNumberOfOutputFeatures();
        int coreSize = lA->getCoreSize();
        bool isConvLayer=true; // are we reading weights for a convolutional or a full connected layer?
        if(inputDepth<1)
        {
            stringstream e;
            e << "Invalid input depth for layer '" << layerName << "' (" << inputDepth << ")";
            THROWERROR(e.str());
        }
        if(outputDepth<1)
        {
            stringstream e;
            e << "Invalid output depth for layer '" << layerName << "' (" << outputDepth << ")";
            THROWERROR(e.str());
        }
        if(coreSize<1)
        {
            isConvLayer=false; // this is a full connected layer
        }
        cout << "### reading weights now..." << endl;


        mpz_class numberOfWeightsRead;
        mpz_class expectedNumberOfWeights;
        if(isConvLayer==true)
        {
            expectedNumberOfWeights = mpz_class(inputDepth)*mpz_class(outputDepth)*mpz_class(coreSize)*mpz_class(coreSize);
            numberOfWeightsRead = readWeightsForConvolutional(file,lA,linebuffer);
        }
        else
        {
            numberOfWeightsRead = readWeightsForFullConnected(file,lA,linebuffer);
            expectedNumberOfWeights = numberOfWeightsRead;
        }

        if(expectedNumberOfWeights!=numberOfWeightsRead)
        {
            stringstream e;
            e << "For Layer '" << layerName << "': Expected number of weights (" << expectedNumberOfWeights
              << ") is not equal to actual number of weights (" << numberOfWeightsRead << ") in file '" << pathToFile << "'";
            THROWERROR(e.str());
        }

        file.close();
    }

    void KerasToNeuralNetwork::readBiasOfCSV(string pathToFile, string layerName)
    {
        ifstream file;
        file.open(pathToFile.c_str());
        if(!file.is_open())
        {
            // pooling for example has no bias file
            return;
        }
        string linebuffer;
        while(!file.eof())
        {
            getline(file,linebuffer);
            if(linebuffer=="")
            {
                continue;
            }
            double bias = stod(linebuffer);
            if(abs(bias) >= 1) cout << "### BIAS >= 1" << endl;
            this->layerArgs[layerName]->addBias(bias);
        }

        file.close();
    }

    void KerasToNeuralNetwork::setUndefinedInputShapes(unsigned int layerPosition)
    {
        if(layerPosition<1)
        {
            THROWERROR("Can't determine input shape of first layer!");
        }
        LayerArguments* lastLA = this->layerArgs[this->orderedLayerNames[layerPosition-1]];
        LayerArguments* actualLA = this->layerArgs[this->orderedLayerNames[layerPosition]];
        if(actualLA->getInputWidth()==-1)
        {
            int inputWidth = 1;
            int inputHeight = 1;
            int inputDepth = lastLA->getNumberOfOutputFeatures();
            if(lastLA->getLayerType()=="Convolutional" || lastLA->getLayerType()=="Pooling")
            {
                /*
                inputWidth = (((lastLA->getInputWidth()+lastLA->getPaddingLeft()+lastLA->getPaddingRight())-lastLA->getCoreSize()+1)/lastLA->getStride());
                inputHeight = (((lastLA->getInputHeight()+lastLA->getPaddingTop()+lastLA->getPaddingBot())-lastLA->getCoreSize()+1)/lastLA->getStride());
                 */
                inputWidth = NeuralNetwork::calculateNewOutputSize(lastLA->getInputWidth(),lastLA->getPaddingLeft(),lastLA->getPaddingRight(),lastLA->getCoreSize(),lastLA->getStride());
                inputHeight = NeuralNetwork::calculateNewOutputSize(lastLA->getInputHeight(),lastLA->getPaddingTop(),lastLA->getPaddingBot(),lastLA->getCoreSize(),lastLA->getStride());
            }
            else if(lastLA->getLayerType()=="FullConnected")
            {
                inputWidth = lastLA->getNumberOfOutputFeatures();
            }

            /*
            if(actualLA->getLayerType()=="FullConnected")
            {
                actualLA->setInputWidth(inputWidth*inputHeight*inputDepth);
            }
            else
            {
                actualLA->setInputWidth(inputWidth);
            }
             */
            actualLA->setInputWidth(inputWidth);
        }
        if(actualLA->getInputHeight()==-1)
        {
            int inputHeight = lastLA->getInputHeight();
            if(lastLA->getLayerType()=="Convolutional" || lastLA->getLayerType()=="Pooling")
            {
                inputHeight = NeuralNetwork::calculateNewOutputSize(lastLA->getInputHeight(),lastLA->getPaddingTop(),lastLA->getPaddingBot(),lastLA->getCoreSize(),lastLA->getStride());
                //inputHeight = (((lastLA->getInputHeight()+lastLA->getPaddingTop()+lastLA->getPaddingBot())-lastLA->getCoreSize())/lastLA->getStride())+1;
            }
            else if(lastLA->getLayerType()=="FullConnected")
                inputHeight = 1;
            /*
            if(actualLA->getLayerType()=="FullConnected")
            {
                actualLA->setInputHeight(1);
            }
            else
            {
                actualLA->setInputHeight(inputHeight);
            }
             */
            actualLA->setInputHeight(inputHeight);
        }
        if(actualLA->getInputDepth()==-1)
        {
            /*
            if(actualLA->getLayerType()=="FullConnected")
            {
                actualLA->setInputDepth(1);
            }
            else
            {
                actualLA->setInputDepth(lastLA->getNumberOfOutputFeatures());
            }
             */
            if(lastLA->getLayerType()=="FullConnected")
                actualLA->setInputDepth(1);
            else
                actualLA->setInputDepth(lastLA->getNumberOfOutputFeatures());
        }

    }

    vector<string> KerasToNeuralNetwork::orderLayerNames(vector<pair<string, string> > edges)
    {
        vector<string> unorderedLayers;
        vector<string> orderedLayers;
        for(auto it : edges)
        {
            string from = it.first;
            string to = it.second;
            bool foundFrom = false;
            bool foundTo = false;
            for(auto it2 : unorderedLayers)
            {
                if(from == it2)
                {
                    foundFrom=true;
                    break;
                }
            }
            for(auto it2 : unorderedLayers)
            {
                if(to == it2)
                {
                    foundTo=true;
                    break;
                }
            }
            if(foundFrom == true && foundTo == true)
            {
                continue;
            }
            if(foundFrom == false)
            {
                unorderedLayers.push_back(from);
            }
            if(foundTo == false)
            {
                unorderedLayers.push_back(to);
            }
        }
        string tempLayerName = "Input";
        orderedLayers.push_back(tempLayerName);
        while(unorderedLayers.size() > orderedLayers.size())
        {
            for(auto it : edges)
            {
                if(it.first == tempLayerName)
                {
                    tempLayerName = it.second;
                    orderedLayers.push_back(tempLayerName);
                    break;
                }
            }
        }
        return orderedLayers;
    }

    mpz_class KerasToNeuralNetwork::readWeightsForFullConnected(ifstream& file, LayerArguments* lA, string firstLineAfterHeader)
    {
        string linebuffer=firstLineAfterHeader;

        /*
        unsigned int neuronCounter = 0;
        unsigned int featureCounter = 0;
        unsigned int rowCounter = 0;
        unsigned int columnCounter = 0;

        const unsigned int maxNeuronCounter = lA->getNumberOfOutputFeatures();
        const unsigned int maxFeatureCounter = lA->getInputDepth();
        const unsigned int maxRowCounter = lA->getInputHeight();
        const unsigned int maxColumnCounter = lA->getInputWidth();
         */

        vector<double> fullConnectedWeights; // temporary weights vector
        // read all the lines and put them into the weights-vector
        while(!file.eof())
        {
            // insert weight into weights-vector
            double weight = stod(linebuffer);
            if(abs(weight)>=1.0) cout << "### WEIGHT >= 1!" << endl;
            fullConnectedWeights.push_back(weight);

            // read the next line for the next iteration
            getline(file,linebuffer);
        }

        lA->setWeights(fullConnectedWeights); // the weights are not ordered at the moment
        mpz_class numberOfWeightsRead=mpz_class(fullConnectedWeights.size());
        cout << "### READ " << numberOfWeightsRead << " WEIGHTS" << endl;
        return numberOfWeightsRead;
    }

    mpz_class KerasToNeuralNetwork::readWeightsForConvolutional(ifstream& file, LayerArguments* lA, string firstLineAfterHeader)
    {
        unsigned int inputDepth = lA->getInputDepth();
        unsigned int outputDepth = lA->getNumberOfOutputFeatures();
        unsigned int coreSize = lA->getCoreSize();
        string linebuffer=firstLineAfterHeader;

        mpz_class numberOfWeightsRead = mpz_class(0);

        unsigned int inputFeatureCounter=0;
        unsigned int outputFeatureCounter=0;
        unsigned int indexCounter=0;
        unsigned int indexIndexCounter=0;

        vector<vector<vector<double>>> ConvWeights; // temporary weights vector
        // read all the lines and put them into the weights-vector
        while(!file.eof())
        {
            // check if we need to alloc more space to the weights-vector
            if(ConvWeights.size() < inputFeatureCounter+1)
            {
                ConvWeights.resize(inputFeatureCounter+1);
            }
            if(ConvWeights[inputFeatureCounter].size() < outputFeatureCounter+1)
            {
                ConvWeights[inputFeatureCounter].resize(outputFeatureCounter+1);
            }
            if(ConvWeights[inputFeatureCounter][outputFeatureCounter].size() < indexCounter+1)
            {
                ConvWeights[inputFeatureCounter][outputFeatureCounter].resize(indexCounter+1);
            }

            // insert weight into weights-vector
            double weight = stod(linebuffer);
            if(abs(weight)>=1.0) cout << "### WEIGHT >= 1!" << endl;
            ConvWeights[inputFeatureCounter][outputFeatureCounter][indexCounter] = weight;

            // read the next line for the next iteration
            numberOfWeightsRead=numberOfWeightsRead+mpz_class(1);
            getline(file,linebuffer);

#if 0 // 1 is the old one
            if(outputFeatureCounter==((unsigned int)outputDepth-1))
            {
                outputFeatureCounter=0;
                if(inputFeatureCounter==((unsigned int)inputDepth-1))
                {
                    inputFeatureCounter=0;

                    // determine index
                    indexCounter = (indexCounter+coreSize);
                    if(indexCounter >= (coreSize*coreSize))
                    {
                        indexIndexCounter++;
                        indexCounter = indexIndexCounter;
                    }
                }
                else
                {
                    inputFeatureCounter++;
                }
            }
            else
            {
                outputFeatureCounter++;
            }
#else
            if(outputFeatureCounter==((unsigned int)outputDepth-1))
            {
                outputFeatureCounter=0;
                if(inputFeatureCounter==((unsigned int)inputDepth-1))
                {
                    inputFeatureCounter=0;

                    // determine index
                    if(indexCounter==(unsigned int)(coreSize*coreSize-1))
                    {
                        indexCounter=0;
                    }
                    else
                    {
                        indexCounter++;
                    }
                }
                else
                {
                    inputFeatureCounter++;
                }
            }
            else
            {
                outputFeatureCounter++;
            }
#endif

        }

        lA->setConvWeights(ConvWeights);
        return numberOfWeightsRead;
    }

    vector<double> KerasToNeuralNetwork::flatten4DTensor(vector<vector<vector<vector<double>>>> t)
    {
        vector<double> returnMe;
        for(auto threeDimIt : t)
        {
            for(auto twoDimIt : threeDimIt)
            {
                for(auto oneDimIt : twoDimIt)
                {
                    for(auto value : oneDimIt)
                    {
                        returnMe.push_back(value);
                    }
                }
            }
        }
        return returnMe;
    }

}//namespace flopoco
#endif
