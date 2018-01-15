// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "PoolingLayer.hpp"

//include other Operators
#include "NeuralNetworks/Utility/WindowShiftRegister.hpp"
#include "NeuralNetworks/Utility/PaddingGenerator.hpp"
#include "PoolingCore.hpp"
#include "NeuralNetworks/ActivationFunctions/ReLU.hpp"


using namespace std;
namespace flopoco {




    PoolingLayer::PoolingLayer(Target* target, int wordSize_, int horizontalSize_, int verticalSize_, int numberOfOutputFeatures_, bool outputFeaturesParallel_, int paddingTop_, string paddingType_, int windowSize_, string activationFunction_, int stride_, int paddingBot_, int paddingLeft_, int paddingRight_) :
        Layer(target) {

        // Operator(target), wordSize(wordSize_), horizontalSize(horizontalSize_), verticalSize(verticalSize_), numberOfOutputFeatures(numberOfOutputFeatures_), outputFeaturesParallel(outputFeaturesParallel_), paddingTop(paddingTop_), paddingType(paddingType_), windowSize(windowSize_), activationFunction(activationFunction_), stride(stride_), paddingBot(paddingBot_), paddingLeft(paddingLeft_), paddingRight(paddingRight_){

        this->myArguments->setWordSize(wordSize_); //
        this->myArguments->setInputWidth(horizontalSize_); //
        this->myArguments->setInputHeight(verticalSize_); //
        this->myArguments->setNumberOfOutputFeatures(numberOfOutputFeatures_); //
        this->myArguments->setInputDepth(numberOfOutputFeatures_); //
        this->myArguments->setOutputFeaturesParallel(outputFeaturesParallel_); //
        this->myArguments->setInputFeaturesParallel(outputFeaturesParallel_); //
        this->myArguments->setActivationFunction(activationFunction_); //
        this->myArguments->setLayerType("Pooling"); //
        this->myArguments->setPaddingTop(paddingTop_); //
        this->myArguments->setPaddingType(paddingType_); //
        this->myArguments->setCoreSize(windowSize_); //
        this->myArguments->setStride(stride_); //
        this->myArguments->setPaddingBot(paddingBot_); //
        this->myArguments->setPaddingLeft(paddingLeft_); //
        this->myArguments->setPaddingRight(paddingRight_); //

        generateVHDLCode(target);
    }

    PoolingLayer::PoolingLayer(Target *target, LayerArguments *args) : Layer(target, args)
    {
        generateVHDLCode(target);
    }

    int PoolingLayer::getNumberOfInstances()
    {
        if(this->myArguments->getOutputFeaturesParallel()==false)
        {
            return 1;
        }
        return this->myArguments->getNumberOfOutputFeatures();
    }

    string PoolingLayer::getOutputSignalName(int feature)
    {
        if(this->myArguments->getNumberOfOutputFeatures()<=feature)
        {
            cout << "PoolingLayer.getOutputSignalName: requested output signal doesn't exist" << endl;
        }
        return "R_"+to_string(feature);
    }

    string PoolingLayer::getInputSignalName(int feature)
    {
        if(this->myArguments->getInputDepth()<=feature)
        {
            cout << "PoolingLayer.getInputSignalName: requested input signal doesn't exist" << endl;
        }
        return "X_"+to_string(feature);
    }

    void PoolingLayer::generateVHDLCode(Target* target)
    {
        this->useNumericStd();

        // throw error
        if(myArguments->getStride()<1)
        {
            stringstream e;
            e << "The stride must be > 0";
            THROWERROR(e);
        }

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="PoolingLayer";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // name of the operator
        ostringstream name;
        name << "PoolingLayer_wordSize_" << myArguments->getWordSize() << "_horizontalSize_" << myArguments->getInputWidth() << "_verticalSize_" << myArguments->getInputHeight() << "_numberOfOutputFeatures_" << myArguments->getNumberOfOutputFeatures() << "_outputFeaturesParallel_" << (myArguments->getOutputFeaturesParallel()==true?"1":"0") << "_padTop_" << myArguments->getPaddingTop() << "_windowSize_" << myArguments->getCoreSize() << "_stride_" << myArguments->getStride() << "_padBot_" << myArguments->getPaddingBot() << "_paddingLeft_" << myArguments->getPaddingLeft() << "_paddingRight_" << myArguments->getPaddingRight() << "_paddingType_" << myArguments->getPaddingType();
        setName(name.str());

        // add inputs and outputs
        addInput("newStep",1);
        addOutput("finished",1);

        int forLimit=(this->myArguments->getOutputFeaturesParallel()==true?this->myArguments->getNumberOfOutputFeatures():1);
        int maxLatency=0;
        int poolingCoreMaxLatencyNumber=0;
        for(int outFC=0; outFC<forLimit; outFC++)
        {

            addInput("X_"+to_string(outFC),myArguments->getWordSize());
            addOutput("R_"+to_string(outFC),myArguments->getWordSize());
            addOutput("getNewData_"+to_string(outFC),1);
            addOutput("validData_o"+to_string(outFC),1);
            addInput("validData_i"+to_string(outFC),1);

            //shift reg
            WindowShiftRegister* winOp = new WindowShiftRegister(target,myArguments->getWordSize(),myArguments->getCoreSize(),myArguments->getInputWidth());
            addSubComponent(winOp);
            inPortMap(winOp,"X","X_"+to_string(outFC));
            inPortMap(winOp,"enable","validData_i"+to_string(outFC));
            inPortMap(winOp,"newStep","newStep");
            for(int i=0; i<myArguments->getCoreSize()*myArguments->getCoreSize(); i++)
            {
                outPortMap(winOp,winOp->outputNames[i],"winShiftOut_"+to_string(i)+"_"+to_string(outFC),true);
            }
            vhdl << instance(winOp,"winShiftInstance_"+to_string(outFC));

            //padding
            PaddingGenerator* padOp = new PaddingGenerator(target,myArguments->getWordSize(),myArguments->getCoreSize(),myArguments->getInputWidth(),myArguments->getInputHeight(),myArguments->getPaddingTop(),myArguments->getStride(),myArguments->getPaddingType(),(outFC==0?true:false),myArguments->getPaddingBot(),myArguments->getPaddingLeft(),myArguments->getPaddingRight());
            addSubComponent(padOp);
            for(int i=0; i<myArguments->getCoreSize()*myArguments->getCoreSize(); i++)
            {
                inPortMap(padOp,"X"+to_string(i),"winShiftOut_"+to_string(i)+"_"+to_string(outFC));
                outPortMap(padOp,"R"+to_string(i),"paddingOut_"+to_string(i)+"_"+to_string(outFC),true);
            }
            inPortMap(padOp,"validData_i","validData_i"+to_string(outFC));
            inPortMap(padOp,"newStep","newStep");
            outPortMap(padOp,"getNewData","getNewData_"+to_string(outFC),false);

            if(outFC==0)
            {
                outPortMap(padOp,"finished","finished",false);
                outPortMap(padOp,"validData_o","validData_temp",true);
            }
            vhdl << instance(padOp,"paddingInstance_"+to_string(outFC));

            //pooling core
            PoolingCore* pooOp = new PoolingCore(target,myArguments->getWordSize(),myArguments->getCoreSize());

            addSubComponent(pooOp);
            for(int i=0; i<myArguments->getCoreSize()*myArguments->getCoreSize(); i++)
            {
                inPortMap(pooOp,"X"+to_string(i),"paddingOut_"+to_string(i)+"_"+to_string(outFC));
            }
            outPortMap(pooOp,"R","poolOut_"+to_string(outFC));
            if(pooOp->getPipelineDepth()>maxLatency)
            {
                poolingCoreMaxLatencyNumber=outFC;
                maxLatency=pooOp->getPipelineDepth();
            }

            vhdl << instance(pooOp,"poolingInstance_"+to_string(outFC));
        }

        if(maxLatency>0)
        {
            syncCycleFromSignal("poolOut_"+to_string(poolingCoreMaxLatencyNumber));
        }

        if(myArguments->getActivationFunction()=="ReLU")
        {
            for(int outFC=0; outFC<forLimit; outFC++)
            {
                ReLU* relOp = new ReLU(target,this->myArguments->getWordSize());
                addSubComponent(relOp);
                inPortMap(relOp,"X","poolOut_"+to_string(outFC));
                outPortMap(relOp,"R","reluOut_"+to_string(outFC),true);
                vhdl << instance(relOp, "reluInstance_"+to_string(outFC));

                vhdl << "R_" << outFC << " <= reluOut_" << outFC << ";" << endl
                     << "validData_o" << outFC << " <= validData_temp;" << endl;
            }
        }
        else
        {
            for(int outFC=0; outFC<forLimit; outFC++)
            {
                vhdl << "R_" << outFC << " <= poolOut_" << outFC << ";" << endl
                     << "validData_o" << outFC << " <= validData_temp;" << endl;
            }
        }
    }

}//namespace flopoco
