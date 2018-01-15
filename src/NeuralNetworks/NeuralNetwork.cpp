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


using namespace std;
namespace flopoco {

    NeuralNetwork::NeuralNetwork(Target* target, map <string, LayerArguments*> layers_, vector <pair <string, string>> edges_, string name_) :
        Operator(target), layerArgs(layers_), edges(edges_), name(name_) {

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="NeuralNetwork";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // definition of the name of the operator
        setName(name);

        // create vhdl code for this network
        buildNeuralNetwork(target);
    }

    NeuralNetwork::NeuralNetwork(Target *target, string name_) : Operator(target), name(name_)
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

    void NeuralNetwork::buildNeuralNetwork(Target* target)
    {
        cout << "###        BEGINNING TO BUILD NETWORK NOW!" << endl;
        Operator* globalOp = buildGlobalController(target);
        cout << "###        GLOBAL CONTROLLER BUILT!" << endl;
        map<string, Layer*> l = buildAllLayers(target,globalOp);
        map<string, Layer*> * lP = &l;
        cout << "###        LAYERS BUILT!" << endl;
        buildAllEdges(target,lP);
        cout << "###        CONNECTIONS BUILT!" << endl;
        instEverything(lP,globalOp);
        cout << "###        VHDL CODE READY!" << endl;
    }

    map<string, Layer*> NeuralNetwork::buildAllLayers(Target* target, Operator* globalOp)
    {
        unsigned int globalOpPortCounter=0;
        map <string, Layer*> returnMe;
        for(auto it : this->layerArgs)
        {
            cout << "###            BUILDING LAYER '" << it.first << "' NOW!" << endl;
            it.second->printLayerArguments();
            if(it.second->getLayerType()!="Input" && it.second->getLayerType()!="Output")
            {
                Layer* op=buildSpecificLayer(target, it.first);
                addSubComponent(op);

                inPortMap(op,"newStep","newStep");
                cout << "### Bop1" << endl;
                string tmp = "finished_"+to_string(globalOpPortCounter);
                declare(tmp,1);
                outPortMap(op,"finished",tmp,false);

                //////////////////////////////////IS THAT SMART?//////////////////////////////////////
                //this->syncCycleFromSignal("finished_"+to_string(globalOpPortCounter));
                //////////////////////////////////////////////////////////////////////////////////////

                cout << "### Bop2" << endl;

                inPortMap(globalOp,tmp,tmp);
                cout << "### Bop3" << endl;
                globalOpPortCounter++;

                returnMe[it.first]=op;
            }
        }
        return returnMe;
    }

    Layer* NeuralNetwork::buildSpecificLayer(Target* target, string nameOfThisLayer)
    {
        cout << "###            IN FUNCTION 'buildSpecificLayer'" << endl;
        LayerArguments* lA = this->layerArgs[nameOfThisLayer];
        cout << "###            GOT LAYER ARGUMENTS" << endl;
        string layerT = lA->getLayerType();
        cout << "###            LAYER TYPE: '" << layerT << "'" << endl;
        Layer* lay=nullptr;
        if(layerT=="Input" || layerT=="Output")
        {
            cout << "###            Input/Output detected!" << endl;
            // do nothing, this will be handled when building edges
        }
        else if(layerT=="Convolutional")
        {
            cout << "###            Convlolutional Layer detected!" << endl;
            lay = buildConvLayer(target, nameOfThisLayer);
        }
        else if(layerT=="Pooling")
        {
            cout << "###            Pooling Layer detected!" << endl;
            lay = buildPoolingLayer(target, nameOfThisLayer);
        }
        else if(layerT=="FullConnected")
        {
            cout << "###            Full Connected Layer detected!" << endl;
            lay = buildFullConnectedLayer(target, nameOfThisLayer);
        }
        else
        {
            stringstream e;
            e << "Undefined Layer Type '" << layerT << "'!";
            cout << e.str() << endl;
            THROWERROR(e.str());
        }
        cout << "###            IN FUNCTION 'buildSpecificLayer', END. RETURNING LAYER NOW" << endl;
        return lay;
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
        this->checkConnection(fromOp, toOp);
        // for each feature, build one memory block
        for(int numIt=0; numIt<fromOp->myArguments->getNumberOfOutputFeatures(); numIt++)
        {
            int addressWidth=ceil(log2(toOp->myArguments->getInputHeight()*toOp->myArguments->getInputWidth())); // that's how much space in BRAM is needed to store one feature
            ControlledMemory* mem = new ControlledMemory(target, fromOp->myArguments->getWordSize(), addressWidth);
            this->addSubComponent(mem);
            this->inPortMap(mem,"newStep","newStep");

            //connect fromLayer with memory block
            string tmp;

            tmp=from+"_out_"+to_string(numIt);
            this->declare(tmp,fromOp->getWidthByPortName(fromOp->getOutputSignalName(numIt)));
            this->outPortMap(fromOp,fromOp->getOutputSignalName(numIt),tmp,false);
            this->inPortMap(mem,"data_i",tmp);

            tmp=from+"_validData_o_"+to_string(numIt);
            this->declare(tmp,fromOp->getWidthByPortName("validData_o_"+to_string(numIt)));
            this->outPortMap(fromOp,"validData_o_"+to_string(numIt),tmp,false);
            this->inPortMap(mem,"validData_i",tmp);




            //connect toLayer with memory block
            tmp=to+"_getNewData_"+to_string(numIt);
            this->declare(tmp,toOp->getWidthByPortName("getNewData_"+to_string(numIt)));
            this->outPortMap(toOp,"getNewData_"+to_string(numIt),tmp,false);
            this->inPortMap(mem,"getNewData",tmp);

            tmp=to+"_validData_i_"+to_string(numIt);
            this->declare(tmp,toOp->getWidthByPortName("validData_i_"+to_string(numIt)));
            this->outPortMap(mem,"validData_o",tmp,false);
            this->inPortMap(toOp,"validData_i_"+to_string(numIt),tmp);

            tmp=to+"_in_"+to_string(numIt);
            this->declare(tmp,toOp->getWidthByPortName(toOp->getInputSignalName(numIt)));
            this->outPortMap(mem,"data_o",tmp,false);
            this->inPortMap(toOp,toOp->getInputSignalName(numIt),tmp);

            this->vhdl << instance(mem,"Memory_from_"+from+"_to_"+to+"_feature_"+to_string(numIt));
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

    void NeuralNetwork::checkConnection(Layer *from, Layer *to)
    {
        // number of features
        if(from->myArguments->getNumberOfOutputFeatures()!=to->myArguments->getInputDepth())
        {
            stringstream e;
            e << "Layer '" << from << "' has " << from->myArguments->getNumberOfOutputFeatures() << " output features and layer '" << to << "' has " << to->myArguments->getInputDepth() << " input features";
            THROWERROR(e.str());
        }

        // data width
        if(from->myArguments->getWordSize()!=to->myArguments->getWordSize())
        {
            stringstream e;
            e << "Layer '" << from << "' has word size " << from->myArguments->getWordSize() << "  and layer '" << to << "' has word size " << to->myArguments->getWordSize();
            THROWERROR(e.str());
        }

        // fraction
        if(from->myArguments->getFraction()!=to->myArguments->getFraction())
        {
            stringstream e;
            e << "Layer '" << from << "' has fraction " << from->myArguments->getFraction() << "  and layer '" << to << "' has fraction " << to->myArguments->getFraction();
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
        cout << "###                GETTING LAYER ARGUMENTS NOW!" << endl;
        LayerArguments* args = this->layerArgs[convLayer];
        cout << "###                GOING INTO CONSTRUCTOR NOW!" << endl;
        ConvolutionalLayer* convL = new ConvolutionalLayer(target,args);
//        ConvolutionalLayer* convL = new ConvolutionalLayer(target,args->getWordSize(),args->getFraction(),args->getInputWidth(),args->getInputHeight(),args->getCoreSize(),args->getNumberOfOutputFeatures(),
//                                                           args->getInputDepth(),args->getConvWeights(),args->getWeightWordSize(),args->getWeightFraction(),args->getPadding(),args->getPaddingType(),false,
//                                                           false,args->getId());
        //this->instMe[convL]=convLayer;
        cout << "###                LAYER '" << convLayer << "' BUILT!" << endl;

        return convL;
    }

    Layer* NeuralNetwork::buildFullConnectedLayer(Target* target, string FCLayer)
    {
        cout << endl << "NEURALNETWORK.BUILDFULLCONNECTEDLAYER: I AM NOT IMPLEMENTED YET!" << endl << endl;
        return nullptr;
    }

    Layer* NeuralNetwork::buildPoolingLayer(Target* target, string poolingLayer)
    {
        LayerArguments* args = this->layerArgs[poolingLayer];
        PoolingLayer* pooL = new PoolingLayer(target,args);
//        PoolingLayer* pooL = new PoolingLayer(target,args->getWordSize(),args->getInputWidth(),args->getInputHeight(),args->getNumberOfOutputFeatures(),
//                                              false,args->getPadding(),args->getPaddingType(),args->getCoreSize());
        //this->instMe[pooL]=poolingLayer;

        return pooL;
    }

    Layer* NeuralNetwork::buildInput(Target* target, string nextLayerName, Layer* nextLayer)
    {
        int inputFeatures = nextLayer->myArguments->getInputDepth();
        int wordSize = nextLayer->myArguments->getWordSize();
        int addressWidth = ceil(log2(nextLayer->myArguments->getInputHeight()*nextLayer->myArguments->getInputWidth()));
        Layer* inputLayer = new InputLayer(target,inputFeatures,wordSize,addressWidth);
        this->addSubComponent(inputLayer);
        this->inPortMap(inputLayer,"newStep","newStep");

        for(int i=0;i<inputFeatures;i++)
        {
            addInput("Data_i_"+to_string(i),wordSize);
            addInput("validData_i_"+to_string(i),1);
            addInput("newDataSet_"+to_string(i),1);

            this->inPortMap(inputLayer,inputLayer->getInputSignalName(i),"Data_i_"+to_string(i));
            this->inPortMap(inputLayer,"validData_i_"+to_string(i),"validData_i_"+to_string(i));
            this->inPortMap(inputLayer,"newDataSet_"+to_string(i),"newDataSet_"+to_string(i));

            string tmp;
            tmp=nextLayerName+"_"+to_string(i);
            this->declare(tmp,nextLayer->getWidthByPortName(nextLayer->getInputSignalName(i)));
            this->outPortMap(inputLayer,inputLayer->getOutputSignalName(i),tmp,false);
            this->inPortMap(nextLayer,nextLayer->getInputSignalName(i),tmp);

            tmp=nextLayerName+"_validData_i_"+to_string(i);
            this->declare(tmp,nextLayer->getWidthByPortName("validData_i_"+to_string(i)));
            this->outPortMap(inputLayer,"validData_o_"+to_string(i),tmp,false);
            this->inPortMap(nextLayer,"validData_i_"+to_string(i),tmp);

            tmp=nextLayerName+"_getNewData_"+to_string(i);
            this->declare(tmp,nextLayer->getWidthByPortName("getNewData_"+to_string(i)));
            this->outPortMap(nextLayer,"getNewData_"+to_string(i),tmp,false);
            this->inPortMap(inputLayer,"getNewData_"+to_string(i),tmp);
        }
        return inputLayer;
    }

    Layer* NeuralNetwork::buildOutput(Target* target, string lastLayerName, Layer* lastLayer)
    {
        cout << "### BUILDING OUTPUT NOW" << endl;
        int howMany = lastLayer->myArguments->getNumberOfOutputFeatures();
        int wordSize = lastLayer->myArguments->getWordSize();
        Layer* outputLayer = new OutputLayer(target,howMany,wordSize);
        addSubComponent(outputLayer);
        this->inPortMap(outputLayer,"newStep","newStep");

        for(int i=0;i<howMany;i++)
        {
            cout << "### i=" << i << endl;
            string tmp;
            tmp=lastLayerName+"_"+to_string(i);
            this->declare(tmp,lastLayer->getWidthByPortName(lastLayer->getOutputSignalName(i)));
            outPortMap(lastLayer,lastLayer->getOutputSignalName(i),tmp,false);
            inPortMap(outputLayer,outputLayer->getInputSignalName(i),tmp);

            tmp="Data_o_"+to_string(i);
            cout << "### Adding Output '" << tmp << "'" << endl;
            //this->declare(tmp,outputLayer->getWidthByPortName(outputLayer->getOutputSignalName(i)));
            addOutput(tmp,outputLayer->getWidthByPortName(outputLayer->getOutputSignalName(i)));
            outPortMap(outputLayer,outputLayer->getOutputSignalName(i),tmp,false);

            tmp=lastLayerName+"_validData_o_"+to_string(i);
            this->declare(tmp,lastLayer->getWidthByPortName("validData_o_"+to_string(i)));
            outPortMap(lastLayer,"validData_o_"+to_string(i),tmp,false);
            inPortMap(outputLayer,"validData_i_"+to_string(i),tmp);
        }
        return outputLayer;
    }

    Operator* NeuralNetwork::buildGlobalController(Target* target)
    {
        unsigned int numberOfLayers = this->layerArgs.size();

        GlobalController* glob = new GlobalController(target,numberOfLayers);
        //this->instMe[glob] = "GlobalController";
        //this->globalOp = glob;

        outPortMap(glob,"newStep","newStep",true);

        return glob;
    }

    void NeuralNetwork::instEverything(map<string, Layer*> *instMe, Operator* globalOp)
    {
        for(auto it : (*instMe))
        {
            cout << "### instantiating '" << it.first << "' now" << endl;
            if(it.second==nullptr) cout << "### FUCK THAT SHIT" << endl;
            this-> vhdl << instance(it.second, it.first+"_instance");
        }

        cout << "### instantiating 'GlobalController' now"  << endl;
        this->vhdl << instance(globalOp, "GlobalController_instance");
    }

}//namespace flopoco
