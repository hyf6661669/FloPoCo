// general c++ library for manipulating streams
#include <iostream>
#include <sstream>
#include <algorithm>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "NeuralNetwork.hpp"

// include all layers (including memory access)
#include "NeuralNetworks/MemoryManagement/ControlledMemory.hpp"
#include "NeuralNetworks/Layers/Input/InputLayer.hpp"
#include "NeuralNetworks/Layers/Output/OutputLayer.hpp"
#include "NeuralNetworks/Layers/Convolutional/ConvolutionalLayer.hpp"
#include "NeuralNetworks/Layers/FullConnected/FullConnectedLayer.hpp"
#include "NeuralNetworks/Layers/Pooling/PoolingLayer.hpp"
#include "NeuralNetworks/Utility/GlobalController.hpp"
#include "NeuralNetworks/MemoryManagement/DmaAccess/GlobalDMAController.hpp"

#define WRITE_WEIGHTS_AS_HEX 0
#define LUT_BASED_ADDR_CALC_LIMIT 10

using namespace std;
namespace flopoco {

    NeuralNetwork::NeuralNetwork(Target* target, map <string, LayerArguments*> layers_, vector <pair <string, string>> edges_, string name_, bool serial, char roundingType_, char multiplicationTypeChar_, bool onlyOutputClassIndex_) :
        Operator(target), layerArgs(layers_), edges(edges_), name(name_), roundingType(roundingType_), multiplicationTypeChar(multiplicationTypeChar_), onlyOutputClassIndex(onlyOutputClassIndex_)
    {

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="NeuralNetwork";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // definition of the name of the operator
        setName(name);

        // create vhdl code for this network
        buildNeuralNetwork(target,serial);
    }

    NeuralNetwork::NeuralNetwork(Target *target, string name_, char roundingType_, char multiplicationTypeChar_, bool onlyOutputClassIndex_) :
            Operator(target), name(name_), roundingType(roundingType_), multiplicationTypeChar(multiplicationTypeChar_), onlyOutputClassIndex(onlyOutputClassIndex_)
    {
        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="NeuralNetwork";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // definition of the name of the operator
        setName(name);

        cout << "Add Layers and Edges. Then you can build the FloPoCo-Operator with NeuralNetwork::buildNeuralNetwork()" << endl;
    }

    void NeuralNetwork::addLayer(string uniqueLayerName, LayerArguments *arguments)
    {
        this->layerArgs[uniqueLayerName]=arguments;
    }

    void NeuralNetwork::addEdge(string from, string to)
    {
        pair<string, string> p;
        p.first=from;
        p.second=to;
        this->edges.push_back(p);
    }

    void NeuralNetwork::setLayers(map<string, LayerArguments *> layers_)
    {
        this->layerArgs = layers_;
    }

    void NeuralNetwork::setEdges(vector<pair<string, string> > edges_)
    {
        this->edges = edges_;
    }

    void NeuralNetwork::buildNeuralNetwork(Target* target, bool serial)
    {
        this->outputsAreValidAlreadyGenerated = false;
        this->fullConnectedWeightsPerDMAAccess = 32;
        checkForErrors();
        handleSerialCalculation(serial);
        Operator* globalOp = buildGlobalController(target);
        map<string, Layer*> l = buildAllLayers(target,globalOp);
        map<string, Layer*> * lP = &l;
        buildAllEdges(target,lP);
        this->checkAndHandleDMA(target,lP);
        instEverything(lP,globalOp);
    }

    void NeuralNetwork::buildOnlyOneLayer(Target *target, string layerName)
    {
        LayerArguments* lA = getLayerArguments(layerName);
        string layerType = lA->getLayerType();
        Layer* layerPtr = nullptr;
        if(layerType=="Input")
        {
            LayerArguments* nextLayer = this->layerArgs[this->findNextLayerName(layerName)];
            int bramAddressWidth = ceil(log2(nextLayer->getInputHeight()*nextLayer->getInputWidth()));
            layerPtr = new InputLayer(target, this, nextLayer->getInputDepth(), nextLayer->getWordSize(), bramAddressWidth, ((!nextLayer->getCalcAllParallel()) && (nextLayer->getLayerType()=="Convolutional")));
        }
        else if(layerType=="Output")
        {
            LayerArguments* lastLayer = this->layerArgs[this->findLastLayerName(layerName)];

            // calculate width
            int leftPad = lastLayer->getPaddingLeft();
            int rightPad = lastLayer->getPaddingRight();
            int topPad = lastLayer->getPaddingTop();
            int botPad = lastLayer->getPaddingBot();
            int stride = lastLayer->getStride();
            int coreS = lastLayer->getCoreSize();
            // check, which padding was set (either left or top must be set)
            if(leftPad<0)
            {
                leftPad=topPad;
            }
            else if(topPad<0)
            {
                topPad=leftPad;
            }
            if(botPad<0)
            {
                botPad=((coreS%2==0)?(topPad-1):topPad);
            }
            if(rightPad<0)
            {
                rightPad=((coreS%2==0)?(leftPad-1):leftPad);
            }

            if((lastLayer->getInputWidth()+leftPad+rightPad-coreS)%stride!=0)
            {
                stringstream e;
                e << "Can't calculate output width of Layer '" << this->findLastLayerName(layerName) << "'";
                THROWERROR(e.str());
            }
            int width = ((lastLayer->getInputWidth()+leftPad+rightPad-coreS)/stride)+1;

            // calculate height
            if((lastLayer->getInputHeight()+topPad+botPad-coreS)%stride!=0)
            {
                stringstream e;
                e << "Can't calculate output height of Layer '" << this->findLastLayerName(layerName) << "'";
                THROWERROR(e.str());
            }
            int height = ((lastLayer->getInputWidth()+topPad+botPad-coreS)/stride)+1;

            layerPtr = new OutputLayer(target, this, lastLayer->getNumberOfOutputFeatures(),
                                       lastLayer->getWordSize(),
                                       true, width, height, true);
        }
        else if(layerType=="Convolutional")
        {
            // arithmetic-based address determination if lut-based address determination wouldn't work because the lut would turn out too large
            if((ceil(log2(lA->getNumberOfOutputFeatures())) + ceil(log2(lA->getInputDepth()))) >= LUT_BASED_ADDR_CALC_LIMIT)
            {
                lA->setLutBasedAddressCalculation(false);
            }
            else
            {
                lA->setLutBasedAddressCalculation(true);
            }
            layerPtr = new ConvolutionalLayer(target, this, lA, true, false, 0x00);
        }
        else if(layerType=="Pooling")
        {
            layerPtr = new PoolingLayer(target, this, lA);
        }
        else if(layerType=="FullConnected")
        {
            // arithmetic-based address determination if lut-based address determination wouldn't work because the lut would turn out too large
            if((ceil(log2(lA->getNumberOfOutputFeatures())) + ceil(log2(lA->getInputDepth()*lA->getInputHeight()*lA->getInputWidth()))) >= LUT_BASED_ADDR_CALC_LIMIT)
            {
                lA->setLutBasedAddressCalculation(false);
            }
            else
            {
                lA->setLutBasedAddressCalculation(true);
            }
            layerPtr = new FullConnectedLayer(target, this, lA, this->fullConnectedWeightsPerDMAAccess, this->roundingType);
        }
        else
        {
            stringstream e;
            e << "Undefined Layer Type '" << layerType << "'!";
            THROWERROR(e.str());
        }
        this->addSubComponent(layerPtr);

        // add inputs and outputs
        vector<Signal*> * ioList = layerPtr->getIOList();
        for(auto signalIt : (*ioList))
        {
            string tmp = signalIt->getName();
            int signalWidth = signalIt->width();
            if(signalIt->type()==Signal::SignalType::in)
            {
                addInput(tmp,signalWidth);
                inPortMap(layerPtr,tmp,tmp);
            }
            else if(signalIt->type()==Signal::SignalType::out)
            {
                addOutput(tmp,signalWidth);
                outPortMap(layerPtr,tmp,tmp,false);
            }
        }

        this->vhdl << instance(layerPtr,layerName+"_instance");
    }

    LayerArguments *NeuralNetwork::getLayerArguments(string layerName)
    {
        for(auto it : this->layerArgs)
        {
            if(it.first==layerName)
            {
                return it.second;
            }
        }
        stringstream e;
        e << "Layer '" << layerName << "' does not exist in this NN!";
        THROWERROR(e.str());
    }

    map<string, Layer*> NeuralNetwork::buildAllLayers(Target* target, Operator* globalOp)
    {
        unsigned int globalOpPortCounter=0;
        map <string, Layer*> returnMe;
        for(auto it : this->layerArgs)
        {
            it.second->printLayerArguments();
            if(it.second->getLayerType()!="Input" && it.second->getLayerType()!="Output")
            {
                Layer* op=buildSpecificLayer(target, it.first);
                addSubComponent(op);

                inPortMap(op,"newStep","newStep");
                string tmp = "finished_"+to_string(globalOpPortCounter);
                declare(tmp,1);
                outPortMap(op,"finished",tmp,false);

                inPortMap(globalOp,tmp,tmp);
                globalOpPortCounter++;

                returnMe[it.first]=op;
            }
        }
        return returnMe;
    }

    Layer* NeuralNetwork::buildSpecificLayer(Target* target, string nameOfThisLayer)
    {
        LayerArguments* lA = this->layerArgs[nameOfThisLayer];
        string layerT = lA->getLayerType();
        Layer* lay=nullptr;
        if(layerT=="Input" || layerT=="Output")
        {
            // do nothing, this will be handled when building edges
        }
        else if(layerT=="Convolutional")
        {
            lay = buildConvLayer(target, nameOfThisLayer);
        }
        else if(layerT=="Pooling")
        {
            lay = buildPoolingLayer(target, nameOfThisLayer);
        }
        else if(layerT=="FullConnected")
        {
            lay = buildFullConnectedLayer(target, nameOfThisLayer);
        }
        else
        {
            stringstream e;
            e << "Undefined Layer Type '" << layerT << "'!";
            THROWERROR(e.str());
        }

        return lay;
    }

    vector<Layer*> NeuralNetwork::needDMA(map<string, Layer *> *layers)
    {
        vector<Layer*> dmaLayers;
        for(map<string, Layer*>::iterator it = layers->begin(); it != layers->end(); ++it)
        {
            LayerArguments* lA = it->second->myArguments;

            // a convolutional layer which is not a parallel implementation
            if(lA->getCalcAllParallel()==false && lA->getLayerType()=="Convolutional")
            {
                dmaLayers.push_back(it->second);
            }
            if(lA->getLayerType()=="FullConnected")
            {
                dmaLayers.push_back(it->second);
            }
        }
        return dmaLayers;
    }

    void NeuralNetwork::checkAndHandleDMA(Target* target, map<string, Layer *> *layers)
    {
        vector<Layer*> dmaLayers = NeuralNetwork::needDMA(layers);
        dmaLayers = this->sortDMALayersByPriority(dmaLayers);
        unsigned int numberOfDMALayers = dmaLayers.size();
        if(numberOfDMALayers > 0)
        {
            // build dma-handler and connect ports to DMA
            GlobalDMAController* dma = new GlobalDMAController(target,numberOfDMALayers);
            addSubComponent(dma);

            // create ports
            // ready to get new weight from DMA
            addOutput(GlobalDMAController::getDMAReadyPortName(),1);
            outPortMap(dma,GlobalDMAController::getDMAReadyPortName(),GlobalDMAController::getDMAReadyPortName(),false);
            // DMA sends valid data
            addInput(GlobalDMAController::getDMAValidPortName(),1);
            inPortMap(dma,GlobalDMAController::getDMAValidPortName(),GlobalDMAController::getDMAValidPortName());
            // DMA sends last data from this packet
            addInput(GlobalDMAController::getDMALastPortName(),1);
            inPortMap(dma,GlobalDMAController::getDMALastPortName(),GlobalDMAController::getDMALastPortName());
            // actual data
            addInput(GlobalDMAController::getDMADataPortName(),32);
            inPortMap(dma,GlobalDMAController::getDMADataPortName(),GlobalDMAController::getDMADataPortName());

            // starts new DMA transfer on rising edge
            addOutput(GlobalDMAController::getNewDMAAccessPortName(),1);
            outPortMap(dma,GlobalDMAController::getNewDMAAccessPortName(),GlobalDMAController::getNewDMAAccessPortName(),false);
            // start address to read data from
            addOutput(GlobalDMAController::getDMAStartAddressPortName(),32);
            outPortMap(dma,GlobalDMAController::getDMAStartAddressPortName(),GlobalDMAController::getDMAStartAddressPortName(),false);
            // number of bytes to read
            addOutput(GlobalDMAController::getDMANumberOfBytesPortName(),23);
            outPortMap(dma,GlobalDMAController::getDMANumberOfBytesPortName(),GlobalDMAController::getDMANumberOfBytesPortName(),false);
            // DMA is programmed for a new read access
            addInput(GlobalDMAController::getDMAFinishedProgrammingPortName(),1);
            inPortMap(dma,GlobalDMAController::getDMAFinishedProgrammingPortName(),GlobalDMAController::getDMAFinishedProgrammingPortName());
            // DMA can be programmed again
            addInput(GlobalDMAController::getDMAReadyToBeProgrammedPortName(),1);
            inPortMap(dma,GlobalDMAController::getDMAReadyToBeProgrammedPortName(),GlobalDMAController::getDMAReadyToBeProgrammedPortName());

            // connect it with the corresponding layers
            for(unsigned int i=0; i<numberOfDMALayers; i++)
            {
                Layer* lay = dmaLayers[i];
                // actual data
                declare(dma->getDataPortName(i)+"_dma", 32);
                outPortMap(dma,dma->getDataPortName(i),dma->getDataPortName(i)+"_dma",false);
                inPortMap(lay,Layer::getWeightInputName(),dma->getDataPortName(i)+"_dma");

                // last data
                declare(dma->getLastDataPortName(i)+"_dma",1);
                outPortMap(dma,dma->getLastDataPortName(i),dma->getLastDataPortName(i)+"_dma",false);
                inPortMap(lay,Layer::getLastWeightsInputName(),dma->getLastDataPortName(i)+"_dma");

                // this layer gets valid data
                declare(Layer::getValidDataName(i)+"_dma",1);
                outPortMap(dma,Layer::getValidDataName(i),Layer::getValidDataName(i)+"_dma",false);
                inPortMap(lay,Layer::getWeightsValidInputName(),Layer::getValidDataName(i)+"_dma");

                // this layer wants to start a new read-access
                declare(dma->getNewReadAccessPortName(i)+"_dma",1);
                inPortMap(dma,dma->getNewReadAccessPortName(i),dma->getNewReadAccessPortName(i)+"_dma");
                outPortMap(lay,Layer::getNewReadAccessOutputName(),dma->getNewReadAccessPortName(i)+"_dma",false);

                // start address this layer wants to read weights from
                declare(dma->getNewStartAddressPortName(i)+"_dma",32);
                inPortMap(dma,dma->getNewStartAddressPortName(i),dma->getNewStartAddressPortName(i)+"_dma");
                outPortMap(lay,Layer::getNextStartAddressOutputName(),dma->getNewStartAddressPortName(i)+"_dma",false);

                // number of bytes this layer wants to read
                declare(dma->getNumberOfTotalBytesPortName(i)+"_dma",23);
                inPortMap(dma,dma->getNumberOfTotalBytesPortName(i),dma->getNumberOfTotalBytesPortName(i)+"_dma");
                outPortMap(lay,Layer::getNumberOfBytesOutputName(),dma->getNumberOfTotalBytesPortName(i)+"_dma",false);

            }

            // inst it
            this->vhdl << instance(dma,"Global_DMA_Controller");
        }
    }

    vector<Layer*> NeuralNetwork::sortDMALayersByPriority(vector<Layer*> layersToSort)
    {
        // sort layers according to "some" metric, that defines its priority for memory access
        map<Layer*, double> scores;
        vector<Layer*> sortedLayers;
        // calculate the metric for each layer and count number of conv layers
        unsigned int numberOfConvLayers=0;
        for(auto it : layersToSort)
        {
            int imageSize = it->myArguments->getInputWidth()*it->myArguments->getInputHeight();
            int numInputFeatures = it->myArguments->getInputDepth();
            int numOutputFeatures = it->myArguments->getNumberOfOutputFeatures();
            int weightsPerMemoryAccess;
            if(it->myArguments->getLayerType()=="Convolutional")
            {
                numberOfConvLayers++;
                weightsPerMemoryAccess = it->myArguments->getCoreSize() * it->myArguments->getCoreSize();
            }
            else
            {
                weightsPerMemoryAccess = this->fullConnectedWeightsPerDMAAccess;
            }
            double scoreTemp = ((double)(imageSize * numInputFeatures * numOutputFeatures))/((double)weightsPerMemoryAccess);
            scores[it] = scoreTemp;
        }

        while(sortedLayers.size() < layersToSort.size())
        {
            double tmpMaxScore = -1.0;
            Layer* tmpLayer = nullptr;
            for(auto it : scores)
            {
                // Convolutional layers have priority over FullConnected layers in terms of memory bandwidth
                if(numberOfConvLayers>0)
                {
                    if(it.second > tmpMaxScore && it.first->myArguments->getLayerType()=="Convolutional")
                    {
                        tmpMaxScore = it.second;
                        tmpLayer = it.first;
                    }
                }
                else
                {
                    if(it.second > tmpMaxScore)
                    {
                        tmpMaxScore = it.second;
                        tmpLayer = it.first;
                    }
                }
            }
            if(tmpLayer->myArguments->getLayerType()=="Convolutional") numberOfConvLayers--;
            sortedLayers.push_back(tmpLayer);
            scores[tmpLayer] = -2.0;
        }
        return sortedLayers;
    }

    void NeuralNetwork::createHeaderForWeightsArray(unsigned int startaddress,unsigned int numberOfDMALayers)
    {
        ofstream file;
        file.open("Weight_Array_Header.h");
        file << "const unsigned int NeuralNetwork_StartAddress = 0x" << hex << startaddress << dec << ";" << endl;
#if WRITE_WEIGHTS_AS_HEX
        file << "const unsigned int NeuralNetwork_Weights [] = {";
#else
        file << "const int NeuralNetwork_Weights [] = {";
#endif
        unsigned int addressCounter = startaddress;
        vector<vector<unsigned int>> weightsForArray;
        for(auto layerIt : this->layerArgs)
        {
            LayerArguments* lA = layerIt.second;
            lA->setStartAddress(addressCounter);

            if(lA->getLayerType() == "Convolutional")
            {
                vector<unsigned int> v = this->getConvolutionalWeightsForHeaderFile(lA,file);
                weightsForArray.push_back(v);
                addressCounter += (v.size()*4);
            }
            if(lA->getLayerType() == "FullConnected")
            {
                vector<unsigned int> v = this->getFullConnectedWeightsForHeaderFile(lA,file);
                weightsForArray.push_back(v);
                addressCounter += (v.size()*4);
            }
        }
        unsigned int numberOfArrayEntries = 0;
        for(unsigned int j=0; j<weightsForArray.size(); j++)
        {
            vector<unsigned int> v = weightsForArray[j];
            for(unsigned int i=0; i<v.size(); i++)
            {
                unsigned int weight = v[i];
#if WRITE_WEIGHTS_AS_HEX
                file << "0x" << hex << weight << dec;
#else
                file << (int)weight;
#endif
                if(j==weightsForArray.size()-1 && i == v.size()-1)
                {
                    file << "};" << endl;
                }
                else
                {
                    file << ", ";
                }
                numberOfArrayEntries++;
            }
        }
        file << "const unsigned int NeuralNetwork_ArraySize = " << numberOfArrayEntries << ";" << endl;

        file << endl;
        file << "void NeuralNetwork_FillDDR()" << endl;
        file << "{" << endl;
        file << tab << "unsigned int i;" << endl;
#if WRITE_WEIGHTS_AS_HEX
        file << tab << "unsigned int* ptr = (unsigned int*) NeuralNetwork_StartAddress;" << endl;
#else
        file << tab << "int* ptr = (int*) NeuralNetwork_StartAddress;" << endl;
#endif
        file << tab << "for(i=0; i<NeuralNetwork_ArraySize; i++)" << endl;
        file << tab << "{" << endl;
        file << tab << tab << "*ptr = NeuralNetwork_Weights[i];" << endl;
        file << tab << tab << "ptr++;" << endl;
        file << tab << "}" << endl;
        file << "}" << endl;
        file.close();
        return;
    }

    vector<unsigned int> NeuralNetwork::getConvolutionalWeightsForHeaderFile(LayerArguments* lA, ofstream &file)
    {
        if(lA->getLayerType() != "Convolutional")
        {
            THROWERROR("Layer Type 'Convolutional' expected, got '"+lA->getLayerType()+"' instead");
        }
        int weightsPerDataBeat = (int)floor(32.0/((double)lA->getWeightWordSize()));
        int vectorIndexCounter = -1;
        if(weightsPerDataBeat<1)
        {
            THROWERROR("Something is wrong with the weight word size (must be in the range [1; 32])");
        }
        vector<unsigned int> weightsForHeader;
        int weightCounter=weightsPerDataBeat;
        for(int outputFeature=0; outputFeature<lA->getNumberOfOutputFeatures(); outputFeature++)
        {
            for(int inputFeature=0; inputFeature<lA->getInputDepth(); inputFeature++)
            {
                for(int index=0; index<lA->getCoreSize()*lA->getCoreSize(); index++)
                {
                    if(weightCounter==weightsPerDataBeat)
                    {
                        weightsForHeader.push_back(0);
                        weightCounter=0;
                        vectorIndexCounter++;
                        if(vectorIndexCounter<0)
                        {
                            THROWERROR("Overflow while creating Weight_Array_Header.h");
                        }
                    }
                    double weight = lA->getConvWeight(inputFeature,outputFeature,index);
                    int scaledWeight = (int)round(weight * pow(2,lA->getWeightFraction()));
                    unsigned int scaledWeight_u = scaledWeight;
                    unsigned int mask = this->getMaskToCutOffBits(lA->getWeightWordSize()) << (weightCounter * lA->getWeightWordSize());
                    unsigned int getMeInVector = ((scaledWeight_u) << (weightCounter * lA->getWeightWordSize())) & mask;
                    weightsForHeader[vectorIndexCounter] = weightsForHeader[vectorIndexCounter] | getMeInVector;
                    weightCounter++;
                }
                weightCounter=weightsPerDataBeat;
            }
        }

        return weightsForHeader;
    }

    vector<unsigned int> NeuralNetwork::getFullConnectedWeightsForHeaderFile(LayerArguments* lA, ofstream &file)
    {
        if(lA->getLayerType() != "FullConnected")
        {
            THROWERROR("Layer Type 'FullConnected' expected, got '"+lA->getLayerType()+"' instead");
        }
        int weightsPerDataBeat = (int)floor(32.0/((double)lA->getWeightWordSize()));
        int vectorIndexCounter = -1;
        if(weightsPerDataBeat<1)
        {
            THROWERROR("Something is wrong with the weight word size (must be in the range [1; 32])");
        }
        vector<unsigned int> weightsForHeader;
        unsigned int numberOfNeurons = lA->getNumberOfOutputFeatures();
        unsigned int weightsPerNeuron = lA->getInputHeight() * lA->getInputWidth() * lA->getInputDepth();
        unsigned int weightsPerMemoryAccess = this->fullConnectedWeightsPerDMAAccess;
        unsigned int weightsPerDataWord = (unsigned int)floor(32.0 / ((double)lA->getWeightWordSize()));
        unsigned int dmaAccessesPerNeuron = (unsigned int)ceil(((double)weightsPerNeuron) / ((double)weightsPerMemoryAccess));
        unsigned int weightsPerLastWeightAccessOfEachNeuron; // the number of weights read in the last reading sequence for each neuron; this variable name is a mess
        if(weightsPerNeuron % weightsPerMemoryAccess == 0) weightsPerLastWeightAccessOfEachNeuron = weightsPerMemoryAccess;
        else weightsPerLastWeightAccessOfEachNeuron = weightsPerNeuron % weightsPerMemoryAccess;

        unsigned int weightsPerLastDataWordOfEachNeuron; // the number of weights that fit into the last data word in the last reading sequence for each neuron; this variable name is a mess as well
        if(weightsPerNeuron % weightsPerDataWord == 0) weightsPerLastDataWordOfEachNeuron = weightsPerDataWord;
        else weightsPerLastDataWordOfEachNeuron = weightsPerNeuron % weightsPerDataWord;

        unsigned int vectorIndex = 0;
        unsigned int weightVectorIndexCounter = -1; // it gets incremented before the first access on the weight vector

        for(unsigned int neuronCounter=0; neuronCounter<numberOfNeurons; neuronCounter++)
        {
            unsigned int inputNeuronCounter = 0;
            for(unsigned int memoryAccessCounter=0; memoryAccessCounter<dmaAccessesPerNeuron; memoryAccessCounter++)
            {

                unsigned int nextForLimit;

                if(memoryAccessCounter==dmaAccessesPerNeuron-1)
                    nextForLimit = weightsPerLastWeightAccessOfEachNeuron;
                else
                    nextForLimit = weightsPerMemoryAccess;

                unsigned int dataIndexCounter = weightsPerDataWord-1;

                for(unsigned int weightCounter = 0; weightCounter < nextForLimit; weightCounter++)
                {
                    if(dataIndexCounter == weightsPerDataWord-1)
                    {
                        dataIndexCounter = 0;
                        weightVectorIndexCounter++;
                        weightsForHeader.push_back(0);
                    }
                    else
                    {
                        dataIndexCounter++;
                    }

                    double weight = lA->getFullConnectedWeight(inputNeuronCounter,neuronCounter);
                    double scaledWeight = weight * pow(2,lA->getWeightFraction());
                    int scaledWeightInt = (int)round(scaledWeight) & NeuralNetwork::getMaskToCutOffBits(lA->getWeightWordSize());

                    unsigned int shiftedWeight = scaledWeightInt << (dataIndexCounter * lA->getWeightWordSize());
                    weightsForHeader[weightVectorIndexCounter] = weightsForHeader[weightVectorIndexCounter] | shiftedWeight;

                    inputNeuronCounter++;
                }

            }
        }


        return weightsForHeader;
    }

    unsigned int NeuralNetwork::getMaskToCutOffBits(unsigned int width)
    {
        if(width >= 32) return -1;
        unsigned int returnMe = 1 << width;
        returnMe--;
        return returnMe;
    }

    void NeuralNetwork::checkForErrors()
    {
        if(this->layerArgs.size()<3)
        {
            stringstream e;
            e << "The Neural Network has less than 3 Layers (needed: Input; >=1 Hidden Layer; Output): " << this->layerArgs.size();
            THROWERROR(e.str());
        }

        bool foundInputLayer = false;
        bool foundOutputLayer = false;
        for(auto it : this->layerArgs)
        {
            if(it.second->getLayerType()=="Input")
            {
                if(foundInputLayer==false)
                {
                    foundInputLayer=true;
                }
                else
                {
                    stringstream e;
                    e << "The Neural Network has more than 1 Input Layer";
                    THROWERROR(e.str());
                }
            }
            if(it.second->getLayerType()=="Output")
            {
                if(foundOutputLayer==false)
                {
                    foundOutputLayer=true;
                }
                else
                {
                    stringstream e;
                    e << "The Neural Network has more than 1 Output Layer";
                    THROWERROR(e.str());
                }
            }
        }
        if(foundInputLayer==false)
        {
            stringstream e;
            e << "The Neural Network has no Input Layer";
            THROWERROR(e.str());
        }
        if(foundOutputLayer==false)
        {
            stringstream e;
            e << "The Neural Network has no Output Layer";
            THROWERROR(e.str());
        }
    }

    void NeuralNetwork::handleSerialCalculation(bool calcSerial)
    {
        if(calcSerial == true)
        {
            unsigned int numberOfDMALayers = 0;
            for(auto it : this->layerArgs)
            {
                if(it.second->getLayerType() == "Convolutional" || it.second->getLayerType() == "FullConnected")
                {
                    numberOfDMALayers++;
                }
                it.second->setCalcAllParallel(false);
            }

            this->createHeaderForWeightsArray(0x0A000000,numberOfDMALayers);
        }
    }

    void NeuralNetwork::buildAllEdges(Target* target, map<string, Layer*> *layerPtrs)
    {
        // check edges for errors
        this->checkEdges();

        for(auto it : this->edges)
        {
            // this only works that way under the assumption that there is no NN, that has input and output directly connected!
            if(this->layerArgs[it.first]->getLayerType()=="Input")
            {
                Layer* in = this->buildInput(target, it.second, (*layerPtrs)[it.second]);
                (*layerPtrs)["Input"]=in;
            }
            else if(this->layerArgs[it.second]->getLayerType()=="Output")
            {
                Layer* out = this->buildOutput(target, it.first, (*layerPtrs)[it.first]);
                (*layerPtrs)["Output"]=out;
            }
            else
            {
                //build a "normal" connection between 2 layers
                buildSpecificEdge(target,it.first,it.second,(*layerPtrs)[it.first],(*layerPtrs)[it.second]);
            }
        }
    }

    void NeuralNetwork::buildSpecificEdge(Target* target, string from, string to, Layer* fromOp, Layer* toOp)
    {
        this->checkConnection(from, to);

        if(fromOp->getOutputMemoryParallelAccess() != toOp->getInputMemoryParallelAccess())
        {
            THROWERROR("Error with parallel memory accesses between layer '"+from+"' and '"+to+"'");
        }
        bool parallelMemoryAccess = fromOp->getOutputMemoryParallelAccess();

        if(parallelMemoryAccess==true)
        {
            // for each feature, build one memory block
            for(int numIt=0; numIt<fromOp->myArguments->getNumberOfOutputFeatures(); numIt++)
            {
                int addressWidth=ceil(log2(toOp->myArguments->getInputHeight()*toOp->myArguments->getInputWidth()));
                bool outputOverwriteValue = (fromOp->myArguments->getLayerType()=="Convolutional" && fromOp->myArguments->getCalcAllParallel()==false);
                bool needGetNewDataReset = (!toOp->myArguments->getCalcAllParallel()) && (toOp->myArguments->getLayerType()=="Convolutional");
                ControlledMemory* mem = new ControlledMemory(target, fromOp->myArguments->getWordSize(), addressWidth, needGetNewDataReset, needGetNewDataReset, false, false);
                this->addSubComponent(mem);
                this->inPortMap(mem,"newStep","newStep");

                /////////////////////////////////////////
                // connect fromLayer with memory block //
                /////////////////////////////////////////
                string tmp;

                tmp=from+"_out_"+to_string(numIt);
                this->declare(tmp,fromOp->getWidthByPortName(fromOp->getOutputSignalName(numIt)));
                this->outPortMap(fromOp,fromOp->getOutputSignalName(numIt),tmp,false);
                this->inPortMap(mem,"data_i",tmp);

                tmp=from+"_"+Layer::getValidDataName(numIt)+"_o";
                this->declare(tmp,fromOp->getWidthByPortName(Layer::getValidDataName(numIt)+"_o"));
                this->outPortMap(fromOp,Layer::getValidDataName(numIt)+"_o",tmp,false);
                this->inPortMap(mem,"validData_i",tmp);

                if(outputOverwriteValue==true)
                {
                    // Data_overwrite signal
                    tmp=from+"_overwrite_"+to_string(numIt);
                    this->declare(tmp,fromOp->getWidthByPortName(fromOp->getOutputSignalName(numIt)));
                    this->outPortMap(mem,"data_overwrite_o",tmp,false);
                    this->inPortMap(fromOp,fromOp->getIntermediateResultName(numIt),tmp);

                    // getNewData_overwrite signal
                    tmp=from+Layer::getGetNewDataName(numIt)+"_overwrite";
                    this->declare(tmp,1);
                    inPortMap(mem,"getNewData_overwrite",tmp);
                    outPortMap(fromOp,Layer::getGetNewDataName(numIt)+"_intermediate",tmp,false);

                    // getNewData_overwrite_reset signal
                    tmp=from+Layer::getGetNewDataName(numIt)+"_overwrite_reset";
                    this->declare(tmp,1);
                    inPortMap(mem,"getNewData_overwrite_reset",tmp);
                    outPortMap(fromOp,Layer::getGetNewDataName(numIt)+"_intermediate_reset",tmp,false);
                }


                ///////////////////////////////////////
                // connect toLayer with memory block //
                ///////////////////////////////////////
                tmp=to+"_"+Layer::getGetNewDataName(numIt);
                this->declare(tmp,toOp->getWidthByPortName(Layer::getGetNewDataName(numIt)));
                this->outPortMap(toOp,Layer::getGetNewDataName(numIt),tmp,false);
                this->inPortMap(mem,"getNewData",tmp);

                if(needGetNewDataReset==true)
                {
                    tmp=to+Layer::getGetNewDataName(numIt)+"_reset";
                    this->declare(tmp,1);
                    this->outPortMap(toOp,Layer::getGetNewDataName(numIt)+"_reset",tmp,false);
                    this->inPortMap(mem,"getNewData_reset",tmp);
                }

                tmp=to+"_"+Layer::getValidDataName(numIt)+"_i";
                this->declare(tmp,toOp->getWidthByPortName(Layer::getValidDataName(numIt)+"_i"));
                this->outPortMap(mem,"validData_o",tmp,false);
                this->inPortMap(toOp,Layer::getValidDataName(numIt)+"_i",tmp);

                tmp=to+"_in_"+to_string(numIt);
                this->declare(tmp,toOp->getWidthByPortName(toOp->getInputSignalName(numIt)));
                this->outPortMap(mem,"data_o",tmp,false);
                this->inPortMap(toOp,toOp->getInputSignalName(numIt),tmp);

                this->vhdl << instance(mem,"Memory_from_"+from+"_to_"+to+"_feature_"+to_string(numIt));
            }
        }
        else
        {
            // build one "big" memory block, where all features are stored successively
            int addressWidth=ceil(log2(toOp->myArguments->getInputHeight()*toOp->myArguments->getInputWidth()*toOp->myArguments->getInputDepth()));
            bool outputOverwriteValue = (fromOp->myArguments->getLayerType()=="Convolutional" && fromOp->myArguments->getCalcAllParallel()==false);
            bool needGetNewDataReset = ((!toOp->myArguments->getCalcAllParallel()) && ((toOp->myArguments->getLayerType()=="Convolutional") || (toOp->myArguments->getLayerType()=="Pooling"))) || (toOp->myArguments->getLayerType() == "FullConnected");
            bool hasSpecificAddressOnGetNewDataReset = (!toOp->myArguments->getCalcAllParallel()) && ((toOp->myArguments->getLayerType()=="Convolutional") || (toOp->myArguments->getLayerType()=="Pooling"));
            bool hasSpecificAddress = (fromOp->myArguments->getLayerType()=="Convolutional" && fromOp->myArguments->getCalcAllParallel()==false);


            ControlledMemory* mem = new ControlledMemory(target, fromOp->myArguments->getWordSize(), addressWidth, outputOverwriteValue, needGetNewDataReset, hasSpecificAddress, hasSpecificAddressOnGetNewDataReset);
            this->addSubComponent(mem);
            this->inPortMap(mem,"newStep","newStep");

            /////////////////////////////////////////
            // connect fromLayer with memory block //
            /////////////////////////////////////////
            string tmp;


            tmp=from+"_out_"+to_string(0);
            this->declare(tmp,fromOp->getWidthByPortName(fromOp->getOutputSignalName(0)));
            this->outPortMap(fromOp,fromOp->getOutputSignalName(0),tmp,false);
            this->inPortMap(mem,"data_i",tmp);

            tmp=from+"_"+Layer::getValidDataName(0)+"_o";
            this->declare(tmp,fromOp->getWidthByPortName(Layer::getValidDataName(0)+"_o"));
            this->outPortMap(fromOp,Layer::getValidDataName(0)+"_o",tmp,false);
            this->inPortMap(mem,"validData_i",tmp);

            if(outputOverwriteValue==true) {
                // Data_overwrite signal
                tmp = from + "_overwrite_" + to_string(0);
                this->declare(tmp, fromOp->getWidthByPortName(fromOp->getOutputSignalName(0)));
                this->outPortMap(mem, "data_overwrite_o", tmp, false);
                this->inPortMap(fromOp, fromOp->getIntermediateResultName(0), tmp);

                // getNewData_overwrite signal
                tmp = from + Layer::getGetNewDataName(0) + "_overwrite";
                this->declare(tmp, 1);
                inPortMap(mem, "getNewData_overwrite", tmp);
                outPortMap(fromOp, Layer::getGetNewDataName(0) + "_intermediate", tmp, false);

                // getNewData_overwrite_reset signal
                tmp = from + Layer::getGetNewDataName(0) + "_overwrite_reset";
                this->declare(tmp, 1);
                inPortMap(mem, "getNewData_overwrite_reset", tmp);
                outPortMap(fromOp, Layer::getGetNewDataName(0) + "_intermediate_reset", tmp, false);


                // reset_address
                tmp = from + "_reset_address";
                this->declare(tmp, fromOp->getAddressWidth());
                inPortMap(mem, "reset_address", tmp);
                outPortMap(fromOp, "reset_address", tmp, false);

            }


            ///////////////////////////////////////
            // connect toLayer with memory block //
            ///////////////////////////////////////
            tmp=to+"_"+Layer::getGetNewDataName(0);
            this->declare(tmp,toOp->getWidthByPortName(Layer::getGetNewDataName(0)));
            this->outPortMap(toOp,Layer::getGetNewDataName(0),tmp,false);
            this->inPortMap(mem,"getNewData",tmp);

            if(needGetNewDataReset==true)
            {
                tmp=to+Layer::getGetNewDataName(0)+"_reset";
                this->declare(tmp,1);
                this->outPortMap(toOp,Layer::getGetNewDataName(0)+"_reset",tmp,false);
                this->inPortMap(mem,"getNewData_reset",tmp);
            }

            if(hasSpecificAddressOnGetNewDataReset==true)
            {
                tmp=to+Layer::getGetNewDataName(0)+"_reset_address";
                this->declare(tmp,addressWidth);
                this->outPortMap(toOp,Layer::getGetNewDataName(0)+"_reset_address",tmp,false);
                this->inPortMap(mem,"getNewData_reset_address",tmp);
            }

            tmp=to+"_"+Layer::getValidDataName(0)+"_i";
            this->declare(tmp,toOp->getWidthByPortName(Layer::getValidDataName(0)+"_i"));
            this->outPortMap(mem,"validData_o",tmp,false);
            this->inPortMap(toOp,Layer::getValidDataName(0)+"_i",tmp);

            tmp=to+"_in_"+to_string(0);
            this->declare(tmp,toOp->getWidthByPortName(toOp->getInputSignalName(0)));
            this->outPortMap(mem,"data_o",tmp,false);
            this->inPortMap(toOp,toOp->getInputSignalName(0),tmp);

            this->vhdl << instance(mem,"Memory_from_"+from+"_to_"+to+"_serialAccess");
        }
    }

    void NeuralNetwork::checkEdges()
    {
        map <string, unsigned int> inputCount;
        for(auto edgeIt : this->edges)
        {
            if(this->checkIfLayerExists(edgeIt.first)==false)
            {
                stringstream e;
                e << "Layer '" << edgeIt.first << "' does not exist!";
                THROWERROR(e.str());
            }
            if(this->checkIfLayerExists(edgeIt.second)==false)
            {
                stringstream e;
                e << "Layer '" << edgeIt.second << "' in the edge-vector does not exist!";
                THROWERROR(e.str());
            }
            ++inputCount[edgeIt.second];

            if(inputCount[edgeIt.second]>1)
            {
                stringstream e;
                e << "Layer '" << edgeIt.second << "' has more than one input, that's not supported yet";
                THROWERROR(e.str());
            }
        }
    }

    void NeuralNetwork::checkConnection(string from, string to)
    {
        if(this->layerArgs[from]->getLayerType()=="FullConnected" && this->layerArgs[to]->getLayerType()=="FullConnected")
        {
            // number of features
            int num1 = this->layerArgs[from]->getNumberOfOutputFeatures();
            int num2 = this->layerArgs[to]->getInputDepth()*this->layerArgs[to]->getInputHeight()*this->layerArgs[to]->getInputWidth();
            if(num1!=num2)
            {
                stringstream e;
                e << "Layer '" << from << "' has " << num1 << " outputs and layer '" << to << "' has " << num2 << "(" << this->layerArgs[to]->getInputDepth() << ", " << this->layerArgs[to]->getInputHeight() <<  ", " << this->layerArgs[to]->getInputWidth() << ") inputs";
                THROWERROR(e.str());
            }
        }
        else if(this->layerArgs[from]->getLayerType()!="FullConnected" && this->layerArgs[to]->getLayerType()!="FullConnected")
        {
            // number of features
            if(this->layerArgs[from]->getNumberOfOutputFeatures()!=this->layerArgs[to]->getInputDepth())
            {
                stringstream e;
                e << "Layer '" << from << "' has " << this->layerArgs[from]->getNumberOfOutputFeatures() << " output features and layer '" << to << "' has " << this->layerArgs[to]->getInputDepth() << " input features";
                THROWERROR(e.str());
            }
        }
        else
        {
            // other checks not implemented yet...
        }


        // data width
        if(this->layerArgs[from]->getWordSize()!=this->layerArgs[to]->getWordSize())
        {
            stringstream e;
            e << "Layer '" << from << "' has word size " << this->layerArgs[from]->getWordSize() << "  and layer '" << to << "' has word size " << this->layerArgs[to]->getWordSize();
            THROWERROR(e.str());
        }

        // fraction
        if(this->layerArgs[from]->getFraction()!=this->layerArgs[to]->getFraction())
        {
            stringstream e;
            e << "Layer '" << from << "' has fraction " << this->layerArgs[from]->getFraction() << "  and layer '" << to << "' has fraction " << this->layerArgs[to]->getFraction();
            THROWERROR(e.str());
        }
    }

    bool NeuralNetwork::checkIfLayerExists(string layerName)
    {
        map<string, LayerArguments *>::iterator first = this->layerArgs.begin();
        map<string, LayerArguments *>::iterator last = this->layerArgs.end();

        while(first!=last)
        {
            if(first->first==layerName) return true;
            ++first;
        }
        return false;
    }

    Layer* NeuralNetwork::buildConvLayer(Target* target, string convLayer)
    {
        LayerArguments* args = this->layerArgs[convLayer];
		// arithmetic-based address determination if lut-based address determination wouldn't work because the lut would turn out too large
		if((ceil(log2(args->getNumberOfOutputFeatures())) + ceil(log2(args->getInputDepth()))) >= LUT_BASED_ADDR_CALC_LIMIT)
		{
			args->setLutBasedAddressCalculation(false);
		}
		else
		{
			args->setLutBasedAddressCalculation(true);
		}
        bool roundAfterConvCore = false;
        bool useBitHeap = true;
        pair <bool, bool> parallelMemoryAccesses = this->getParallelMemoryAccesses(convLayer);
        ConvolutionalLayer* convL = new ConvolutionalLayer(target,this,args,useBitHeap,roundAfterConvCore,this->roundingType,parallelMemoryAccesses.first,parallelMemoryAccesses.second);
        return convL;
    }

    Layer* NeuralNetwork::buildFullConnectedLayer(Target* target, string FCLayer)
    {
        LayerArguments* args = this->layerArgs[FCLayer];
        // arithmetic-based address determination if lut-based address determination wouldn't work because the lut would turn out too large
		if((ceil(log2(args->getNumberOfOutputFeatures())) + ceil(log2(args->getInputDepth()))) >= LUT_BASED_ADDR_CALC_LIMIT)
		{
			args->setLutBasedAddressCalculation(false);
		}
		else
		{
			args->setLutBasedAddressCalculation(true);
		}
        FullConnectedLayer* fc = new FullConnectedLayer(target,this,args,this->fullConnectedWeightsPerDMAAccess,this->roundingType);
        return fc;
    }

    Layer* NeuralNetwork::buildPoolingLayer(Target* target, string poolingLayer)
    {
        pair <bool, bool> parallelMemoryAccesses = this->getParallelMemoryAccesses(poolingLayer);
        LayerArguments* args = this->layerArgs[poolingLayer];
        PoolingLayer* pooL = new PoolingLayer(target,this,args,parallelMemoryAccesses.first,parallelMemoryAccesses.second);
        return pooL;
    }

    Layer* NeuralNetwork::buildInput(Target* target, string nextLayerName, Layer* nextLayer)
    {
        nextLayer->myArguments->printLayerArguments();
        bool needGetNewDataReset = ((!nextLayer->myArguments->getCalcAllParallel()) && (nextLayer->myArguments->getLayerType()=="Convolutional")) || (nextLayer->myArguments->getLayerType()=="FullConnected");
        int inputFeatures = nextLayer->myArguments->getInputDepth();
        int wordSize = nextLayer->myArguments->getWordSize();
        int addressWidth = ceil(log2(nextLayer->myArguments->getInputHeight()*nextLayer->myArguments->getInputWidth()));
        Layer* inputLayer = new InputLayer(target,this,inputFeatures,wordSize,addressWidth,needGetNewDataReset);
        this->addSubComponent(inputLayer);
        this->inPortMap(inputLayer,"newStep","newStep");

        for(int i=0;i<inputFeatures;i++)
        {
            addInput("Data_i_"+to_string(i),wordSize);
            addInput(Layer::getValidDataName(i)+"_i",1);
            addInput("newDataSet_"+to_string(i),1);

            this->inPortMap(inputLayer,inputLayer->getInputSignalName(i),"Data_i_"+to_string(i));
            this->inPortMap(inputLayer,Layer::getValidDataName(i)+"_i",Layer::getValidDataName(i)+"_i");
            this->inPortMap(inputLayer,"newDataSet_"+to_string(i),"newDataSet_"+to_string(i));

            string tmp;
            tmp=nextLayerName+"_in_"+to_string(i);
            this->declare(tmp,nextLayer->getWidthByPortName(nextLayer->getInputSignalName(i)));
            this->outPortMap(inputLayer,inputLayer->getOutputSignalName(i),tmp,false);
            this->inPortMap(nextLayer,nextLayer->getInputSignalName(i),tmp);

            tmp=nextLayerName+"_"+Layer::getValidDataName(i)+"_i";
            this->declare(tmp,nextLayer->getWidthByPortName(Layer::getValidDataName(i)+"_i"));
            this->outPortMap(inputLayer,Layer::getValidDataName(i)+"_o",tmp,false);
            this->inPortMap(nextLayer,Layer::getValidDataName(i)+"_i",tmp);

            tmp=nextLayerName+"_"+Layer::getGetNewDataName(i);
            this->declare(tmp,nextLayer->getWidthByPortName(Layer::getGetNewDataName(i)));
            this->outPortMap(nextLayer,Layer::getGetNewDataName(i),tmp,false);
            this->inPortMap(inputLayer,Layer::getGetNewDataName(i),tmp);

            if(needGetNewDataReset==true)
            {
                tmp=nextLayerName+"_"+Layer::getGetNewDataName(i)+"_reset";
                this->declare(tmp,1);
                this->outPortMap(nextLayer,Layer::getGetNewDataName(i)+"_reset",tmp,false);
                this->inPortMap(inputLayer,Layer::getGetNewDataName(i)+"_reset",tmp);
            }
        }
        return inputLayer;
    }

    Layer* NeuralNetwork::buildOutput(Target* target, string lastLayerName, Layer* lastLayer)
    {
        OutputLayer* outputLayer = nullptr;

        int howMany = lastLayer->myArguments->getNumberOfOutputFeatures();
        int wordSize = lastLayer->myArguments->getWordSize();

        if(lastLayer->myArguments->getLayerType()=="FullConnected")
        {
            outputLayer = new OutputLayer(target,this,howMany,wordSize,0,0,0,0,0,this->onlyOutputClassIndex);
            addSubComponent(outputLayer);

            /////////////////////////////
            // connect with last layer //
            /////////////////////////////
            string tmp;

            // data
            tmp=lastLayerName+"_out";
            this->declare(tmp,lastLayer->getWidthByPortName(lastLayer->getOutputSignalName(0)));
            outPortMap(lastLayer,lastLayer->getOutputSignalName(0),tmp,false);
            inPortMap(outputLayer,outputLayer->getInputSignalName(0),tmp);

            // validData
            tmp=lastLayerName+"validData_out";
            this->declare(tmp,1);
            outPortMap(lastLayer,lastLayer->getValidDataName(0)+"_o",tmp,false);
            inPortMap(outputLayer,outputLayer->getValidDataName(0)+"_i",tmp);

            ///////////////////////////////////
            // add outputs to neural network //
            ///////////////////////////////////
            if(this->onlyOutputClassIndex==false)
            {
                for(int i=0; i<howMany; i++)
                {
                    tmp="R_"+to_string(i);
                    this->addOutput(tmp,wordSize);
                    outPortMap(outputLayer,outputLayer->getOutputSignalName(i),tmp,false);
                }
            }

            tmp="ClassIndex";
            this->addOutput(tmp,outputLayer->getIndexWordSize());
            outPortMap(outputLayer,"I_max",tmp,false);

            /////////////////////////////////
            // set flag: outputs_are_valid //
            /////////////////////////////////
            outPortMap(outputLayer,"Outputs_Are_Valid","Outputs_Are_Valid",false);
            this->outputsAreValidAlreadyGenerated=true;
        }
        else
        {
            // the output values are valid, when all layers finished their calculation (= newStep is assigned by global controller)
            vhdl << "Outputs_Are_Valid <= newStep;" << endl;
            // calculate width
            int topPad = lastLayer->myArguments->getPaddingTop();
            int botPad = lastLayer->myArguments->getPaddingBot();
            int leftPad = lastLayer->myArguments->getPaddingLeft();
            int rightPad = lastLayer->myArguments->getPaddingRight();
            int stride = lastLayer->myArguments->getStride();
            int coreS = lastLayer->myArguments->getCoreSize();
            // check, which padding was set! (either left or top pad must be set)
            if(leftPad<0)
            {
                leftPad=topPad;
            }
            else if(topPad<0)
            {
                topPad=leftPad;
            }
            if(botPad<0)
            {
                botPad=((coreS%2==0)?(topPad-1):topPad);
            }
            if(rightPad<0)
            {
                rightPad=((coreS%2==0)?(leftPad-1):leftPad);
            }

            int width = NeuralNetwork::calculateNewOutputSize(lastLayer->myArguments->getInputWidth(),leftPad,rightPad,coreS,stride);
            int height = NeuralNetwork::calculateNewOutputSize(lastLayer->myArguments->getInputHeight(),topPad,botPad,coreS,stride);
            bool needIntermediateValue = (!lastLayer->myArguments->getCalcAllParallel()) && lastLayer->myArguments->getLayerType()=="Convolutional";
            outputLayer = new OutputLayer(target,this,howMany,wordSize,true,width,height,true,needIntermediateValue);
            this->addSubComponent(outputLayer);
            for(int i=0;i<howMany;i++)
            {
                string tmp;
                tmp=lastLayerName+"_out_"+to_string(i);
                this->declare(tmp,lastLayer->getWidthByPortName(lastLayer->getOutputSignalName(i)));
                outPortMap(lastLayer,lastLayer->getOutputSignalName(i),tmp,false);
                inPortMap(outputLayer,outputLayer->getInputSignalName(i),tmp);

                tmp="Data_o_"+to_string(i);
                addOutput(tmp,outputLayer->getWidthByPortName(outputLayer->getOutputSignalName(i)));
                outPortMap(outputLayer,outputLayer->getOutputSignalName(i),tmp,false);

                tmp=lastLayerName+"_"+Layer::getValidDataName(i)+"_o";
                this->declare(tmp,lastLayer->getWidthByPortName(Layer::getValidDataName(i)+"_o"));
                outPortMap(lastLayer,Layer::getValidDataName(i)+"_o",tmp,false);
                inPortMap(outputLayer,Layer::getValidDataName(i)+"_i",tmp);

                tmp=Layer::getValidDataName(i)+"_o";
                addOutput(tmp,1);
                outPortMap(outputLayer,tmp,tmp,false);

                tmp="getNewDataFromOutport_"+to_string(i);
                addInput(tmp,1);
                inPortMap(outputLayer,Layer::getGetNewDataName(i),tmp);

                tmp="getNewDataSetFromOutport_"+to_string(i);
                addInput(tmp,1);
                inPortMap(outputLayer,"getNewDataSet_"+to_string(i),tmp);

                if(needIntermediateValue==true)
                {
                    tmp="Intermediate_Result_From_Output_"+to_string(i);
                    outPortMap(outputLayer,outputLayer->getIntermediateResultName(i),tmp,true);
                    inPortMap(lastLayer,lastLayer->getIntermediateResultName(i),tmp);

                    tmp="Get_New_Data_intermediate_From_Output_"+to_string(i);
                    outPortMap(lastLayer,lastLayer->getGetNewDataName(i)+"_intermediate",tmp);
                    inPortMap(outputLayer,outputLayer->getGetNewDataName(i)+"_intermediate",tmp);

                    tmp="Reset_Get_New_Data_intermediate_From_Output_"+to_string(i);
                    outPortMap(lastLayer,Layer::getGetNewDataName(i)+"_intermediate_reset",tmp);
                    inPortMap(outputLayer,outputLayer->getGetNewDataName(i)+"_intermediate_reset",tmp);
                }
            }
        }

        this->inPortMap(outputLayer,"newStep","newStep");



        return outputLayer;
    }

    Operator* NeuralNetwork::buildGlobalController(Target* target)
    {
        unsigned int numberOfFinishedSignals = this->layerArgs.size()-2; // input and output layer don't need a finishedSignal
        GlobalController* glob = new GlobalController(target,numberOfFinishedSignals);
        addSubComponent(glob);
        outPortMap(glob,"newStep","newStep",true);
        this->addOutput("Outputs_Are_Valid",1);
        return glob;
    }

    string NeuralNetwork::findLastLayerName(string actualLayerName) const
    {
        for(auto it : this->edges)
        {
            if(it.second == actualLayerName)
            {
                return it.first;
            }
        }
        stringstream e;
        e << "Layer '" << actualLayerName << "' does not have any preceding layer";
        THROWERROR(e.str());
    }

    string NeuralNetwork::findNextLayerName(string actualLayerName) const
    {
        for(auto it : this->edges)
        {
            if(it.first == actualLayerName)
            {
                return it.second;
            }
        }
        stringstream e;
        e << "Layer '" << actualLayerName << "' does not have any next layer";
        THROWERROR(e.str());
    }

    pair<bool, bool> NeuralNetwork::getParallelMemoryAccesses(string layerName) const
    {
        bool inputMemoryParallelAccess=true;
        bool outputMemoryParallelAccess=true;
        LayerArguments* lALastLayer = this->layerArgs.at(findLastLayerName(layerName));
        LayerArguments* lANextLayer = this->layerArgs.at(findNextLayerName(layerName));
        if(((lALastLayer->getLayerType()=="Convolutional" || lALastLayer->getLayerType()=="Pooling") && lALastLayer->getCalcAllParallel()==false) || (lALastLayer->getLayerType()=="FullConnected"))
        {
            inputMemoryParallelAccess=false;
        }
        if(((lANextLayer->getLayerType()=="Convolutional" || lANextLayer->getLayerType()=="Pooling") && lANextLayer->getCalcAllParallel()==false) || (lANextLayer->getLayerType()=="FullConnected"))
        {
            outputMemoryParallelAccess=false;
        }
        return {inputMemoryParallelAccess, outputMemoryParallelAccess};
    }

    void NeuralNetwork::instEverything(map<string, Layer*> *instMe, Operator* globalOp)
    {
        for(auto it : (*instMe))
        {
            this-> vhdl << instance(it.second, it.first+"_instance");
        }

        this->vhdl << instance(globalOp, "GlobalController_instance");
    }

    unsigned int NeuralNetwork::calculateNewOutputSize(unsigned int oldSize, unsigned int padding1, unsigned int padding2, unsigned int coreSize, unsigned int stride)
    {
        if(((oldSize+padding1+padding2)-coreSize)%stride!=0)
        {
            cout << "Parameters might be wrong, stride (" << stride << ") doesn't fit the values for image size (" << oldSize
                 << "), padding (" << padding1 << ", " << padding2 << "), core size (" << coreSize << ")" << endl;
        }
        unsigned int value = (unsigned int)floor(((double)(oldSize+padding1+padding2-coreSize))/((double)stride))+1;
        return value;
    }

    multiplicationType NeuralNetwork::getMultiplicationType() const
    {
        switch(this->multiplicationTypeChar)
        {
            case 0x00: return multiplicationType::simple; break;
            default: return multiplicationType::none; break;
        }
    }

    string NeuralNetwork::convertMpzToString(mpz_class number, unsigned int width)
    {
        stringstream returnMe;
        unsigned int howMuchShift;
        for(unsigned int i=0; i<width; i++)
        {
            howMuchShift = width-i-1;
            returnMe << ((number >> howMuchShift) & 1);
        }
        return returnMe.str();
    }


}//namespace flopoco
