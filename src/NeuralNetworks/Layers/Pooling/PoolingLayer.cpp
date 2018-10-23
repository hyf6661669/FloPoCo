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
#include "IntAddSubCmp/IntComparator.hpp"
#include "NeuralNetworks/Utility/ModuloCounter.hpp"
#include "PrimitiveComponents/GenericMux.hpp"
#include "PrimitiveComponents/GenericLut.hpp"
#include "NeuralNetworks/Utility/DataGuard.hpp"
#include "NeuralNetworks/Utility/Register.hpp"


using namespace std;
namespace flopoco {




    PoolingLayer::PoolingLayer(Target* target, NeuralNetwork* parent, int wordSize_, int horizontalSize_, int verticalSize_, int numberOfOutputFeatures_, bool calcAllParallel_, int paddingTop_, string paddingType_, int windowSize_, string activationFunction_, int stride_, int paddingBot_, int paddingLeft_, int paddingRight_, bool inputMemoryParallelAccess_, bool outputMemoryParallelAccess_) :
        Layer(target,parent){

        this->myArguments->setWordSize(wordSize_); //
        this->myArguments->setInputWidth(horizontalSize_); //
        this->myArguments->setInputHeight(verticalSize_); //
        this->myArguments->setNumberOfOutputFeatures(numberOfOutputFeatures_); //
        this->myArguments->setInputDepth(numberOfOutputFeatures_); //
        this->myArguments->setCalcAllParallel(calcAllParallel_); //
        this->myArguments->setActivationFunction(activationFunction_); //
        this->myArguments->setLayerType("Pooling"); //
        this->myArguments->setPaddingTop(paddingTop_); //
        this->myArguments->setPaddingType(paddingType_); //
        this->myArguments->setCoreSize(windowSize_); //
        this->myArguments->setStride(stride_); //
        this->myArguments->setPaddingBot(paddingBot_); //
        this->myArguments->setPaddingLeft(paddingLeft_); //
        this->myArguments->setPaddingRight(paddingRight_); //

        this->inputMemoryParallelAccess = inputMemoryParallelAccess_;
        this->outputMemoryParallelAccess = outputMemoryParallelAccess_;

        generateVHDLCode(target);
    }

    PoolingLayer::PoolingLayer(Target *target, NeuralNetwork* parent, LayerArguments *args, bool inputMemoryParallelAccess_, bool outputMemoryParallelAccess_) :
            Layer(target, args, parent)
    {
        this->inputMemoryParallelAccess = inputMemoryParallelAccess_;
        this->outputMemoryParallelAccess = outputMemoryParallelAccess_;

        generateVHDLCode(target);
    }

    int PoolingLayer::getNumberOfInstances()
    {
        if(this->myArguments->getCalcAllParallel()==false)
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
            THROWERROR(e.str());
        }

        if(this->myArguments->getPaddingBot()<0) // padBot
        {
            this->myArguments->setPaddingBot((this->myArguments->getCoreSize()%2==1 || this->myArguments->getPadding()==0)?this->myArguments->getPadding():this->myArguments->getPadding()-1); // padBot=((windowSize%2)==1?padTop:padTop-1);
        }
        if(this->myArguments->getPaddingLeft()<0) // padLeft
        {
            this->myArguments->setPaddingLeft(this->myArguments->getPadding()); // padLeft=padTop;
        }
        if(this->myArguments->getPaddingRight()<0) // padRight
        {
            this->myArguments->setPaddingRight(((this->myArguments->getCoreSize()%2)==1 || this->myArguments->getPadding()==0)?this->myArguments->getPadding():this->myArguments->getPadding()-1); // padRight=((windowSize%2)==1?padTop:padTop-1);
        }

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="PoolingLayer";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // name of the operator
        ostringstream name;
        name << "PoolingLayer_wordSize_" << myArguments->getWordSize() << "_horizontalSize_" << myArguments->getInputWidth() << "_verticalSize_" << myArguments->getInputHeight() << "_numberOfOutputFeatures_" << myArguments->getNumberOfOutputFeatures() << "_outputFeaturesParallel_" << (myArguments->getCalcAllParallel()==true?"1":"0") << "_padTop_" << myArguments->getPaddingTop() << "_windowSize_" << myArguments->getCoreSize() << "_stride_" << myArguments->getStride() << "_padBot_" << myArguments->getPaddingBot() << "_paddingLeft_" << myArguments->getPaddingLeft() << "_paddingRight_" << myArguments->getPaddingRight() << "_paddingType_" << myArguments->getPaddingType();
        setName(name.str());

        if(this->myArguments->getCalcAllParallel()==true && (this->inputMemoryParallelAccess==false || this->outputMemoryParallelAccess==false))
        {
            THROWERROR("This Layer needs parallel access to input/output feature BRAMs");
        }

        // add inputs and outputs
        addInput("newStep",1);
        addOutput("finished",1);

        int forLimit=this->getNumberOfInstances();
        int maxLatency=0;
        int poolingCoreMaxLatencyNumber=0;
        for(int outFC=0; outFC<forLimit; outFC++)
        {
            if(this->myArguments->getCalcAllParallel()==true)
            {
                addInput("X_"+to_string(outFC),myArguments->getWordSize());
                addOutput("R_"+to_string(outFC),myArguments->getWordSize());
                addOutput(Layer::getValidDataName(outFC)+"_o",1);
                addInput(Layer::getValidDataName(outFC)+"_i",1);
                addOutput(Layer::getGetNewDataName(outFC),1);
            }
            else
            {
                // Ports
                if(this->inputMemoryParallelAccess==true)
                {
                    for(int i=0; i<this->myArguments->getNumberOfOutputFeatures(); i++)
                    {
                        addInput("X_"+to_string(i),myArguments->getWordSize());
                        addInput(Layer::getValidDataName(i)+"_i",1);
                        addOutput(Layer::getGetNewDataName(i),1);
                    }
                }
                else
                {
                    addInput("X_"+to_string(0),myArguments->getWordSize());
                    addInput(Layer::getValidDataName(0)+"_i",1);
                    addOutput(Layer::getGetNewDataName(0),1);
                }
                if(this->outputMemoryParallelAccess==true)
                {
                    for(int i=0; i<this->myArguments->getNumberOfOutputFeatures(); i++)
                    {
                        addOutput("R_"+to_string(i),myArguments->getWordSize());
                        addOutput(Layer::getValidDataName(i)+"_o",1);
                    }
                }
                else
                {
                    addOutput("R_"+to_string(0),myArguments->getWordSize());
                    addOutput(Layer::getValidDataName(0)+"_o",1);
                }

                // Feature-Counter
                declare("finished_padding",1);
                declare("finished_padding_reg",1);
                vhdl << declare("finished_padding_rise",1) << " <= finished_padding and (not finished_padding_reg);" << endl;
                vhdl << declare("counter_enable",1) << " <= newStep or (finished_padding_rise and not featureCounter_max);" << endl;
                Register* reg = new Register(target,1,1);
                addSubComponent(reg);
                inPortMap(reg,"X","finished_padding");
                outPortMap(reg,"R","finished_padding_reg",false);
                vhdl << instance(reg,"Finished_Padding_Reg_Instance");

                if(this->myArguments->getNumberOfOutputFeatures()>1)
                {
                    ModuloCounter* modC = new ModuloCounter(target, this->myArguments->getNumberOfOutputFeatures());
                    this->addSubComponent(modC);
                    this->inPortMap(modC,"manualReset","newStep");
                    this->outPortMap(modC,"counter","featureCounter");
                    this->inPortMap(modC,"enable","counter_enable");
                    this->vhdl << instance(modC,"FeatureCounter_instance");
                }
                else
                {
                    this->vhdl << declare("featureCounter",1) << " <= \"0\";" << endl;
                }

                if(this->inputMemoryParallelAccess==false)
                {
                    // create a reset for the input memory since the layer is reading a bit more than it processes every time (that makes padding MUCH easier)
                    this->createInputMemoryReset(target);

                }

            }

            //shift reg
            WindowShiftRegister* winOp = new WindowShiftRegister(target,myArguments->getWordSize(),myArguments->getCoreSize(),myArguments->getInputWidth());
            addSubComponent(winOp);

            if(this->myArguments->getCalcAllParallel()==true)
            {
                inPortMap(winOp,"X","X_"+to_string(outFC));
                inPortMap(winOp,"enable",Layer::getValidDataName(outFC)+"_i");
                inPortMap(winOp,"newStep","newStep");
            }
            else
            {
                if(this->inputMemoryParallelAccess==true)
                {
                    // Mux for all Inputs (X)
                    GenericMux* inputMux = new GenericMux(target, this->myArguments->getWordSize(), this->myArguments->getInputDepth());
                    this->addSubComponent(inputMux);
                    for(unsigned int i=0;i<(unsigned int)this->myArguments->getInputDepth();i++)
                    {
                        this->inPortMap(inputMux,inputMux->getInputName(i),this->getInputSignalName(i));
                    }
                    this->inPortMap(inputMux,inputMux->getSelectName(),"featureCounter");
                    this->outPortMap(inputMux,inputMux->getOutputName(),"X_Mux",true);
                    this->vhdl << instance(inputMux,"InputMux_instance");
                    this->inPortMap(winOp,"X","X_Mux");

                    // Mux for all Valid-signals
                    GenericMux* inputMux2 = new GenericMux(target, 1, this->myArguments->getInputDepth());
                    this->addSubComponent(inputMux2);
                    for(unsigned int i=0;i<(unsigned int)this->myArguments->getInputDepth();i++)
                    {
                        this->inPortMap(inputMux2,inputMux2->getInputName(i),Layer::getValidDataName(i)+"_i");
                    }
                    this->inPortMap(inputMux2,inputMux2->getSelectName(),"featureCounter");
                    this->outPortMap(inputMux2,inputMux2->getOutputName(),"validData_i_Mux",true);
                    this->vhdl << instance(inputMux2,"ValidData_i_Mux_instance");
                    this->inPortMap(winOp,"enable","validData_i_Mux");
                }
                else
                {
                    this->inPortMap(winOp,"X","X_"+to_string(0));
                    this->inPortMap(winOp,"enable",Layer::getValidDataName(0)+"_i");
                }

                // newStep
                this->declare("newStepOrFinishedFeature",1);
                this->inPortMap(winOp,"newStep","newStepOrFinishedFeature");
            }

            for(int i=0; i<myArguments->getCoreSize()*myArguments->getCoreSize(); i++)
            {
                outPortMap(winOp,winOp->outputNames[i],"winShiftOut_"+to_string(i)+"_"+to_string(outFC),true);
            }
            vhdl << instance(winOp,"winShiftInstance_"+to_string(outFC));

            //padding
            bool buildForSerialCalc = false;
            PaddingGenerator* padOp = new PaddingGenerator(target,myArguments->getWordSize(),myArguments->getCoreSize(),myArguments->getInputWidth(),myArguments->getInputHeight(),myArguments->getPaddingTop(),myArguments->getStride(),myArguments->getPaddingType(),(outFC==0?true:false),buildForSerialCalc,myArguments->getPaddingBot(),myArguments->getPaddingLeft(),myArguments->getPaddingRight());
            addSubComponent(padOp);
            for(int i=0; i<myArguments->getCoreSize()*myArguments->getCoreSize(); i++)
            {
                inPortMap(padOp,"X"+to_string(i),"winShiftOut_"+to_string(i)+"_"+to_string(outFC));
                outPortMap(padOp,"R"+to_string(i),"paddingOut_"+to_string(i)+"_"+to_string(outFC),true);
            }

            if(this->myArguments->getCalcAllParallel()==true)
            {
                inPortMap(padOp,"validData_i",Layer::getValidDataName(outFC)+"_i");
                inPortMap(padOp,"newStep","newStep");
                outPortMap(padOp,"getNewData",Layer::getGetNewDataName(outFC),false);
                if(outFC==0)
                {
                    outPortMap(padOp,"finished","finished",false);
                    outPortMap(padOp,"validData_o","validData_temp",true);
                }
            }
            else
            {
                if(this->inputMemoryParallelAccess==true)
                {
                    // validData_i
                    inPortMap(padOp,"validData_i","validData_i_Mux");

                    // Guard for getNewData
                    outPortMap(padOp,"getNewData","getNewData_PadGen");
                    for(unsigned int i=0;i<(unsigned int)this->myArguments->getInputDepth();i++)
                    {
                        DataGuard* datG = new DataGuard(target, 1, this->getSignalByName("featureCounter")->width(),i);
                        this->addSubComponent(datG);
                        this->inPortMap(datG,"Guard_in","featureCounter");
                        this->inPortMap(datG,"Data_in","getNewData_PadGen");
                        this->outPortMap(datG,"Data_out",Layer::getGetNewDataName(i),false);
                        this->vhdl << instance(datG,"GetNewDataGuard_"+to_string(i));
                    }
                }
                else
                {
                    inPortMap(padOp,"validData_i",Layer::getValidDataName(0)+"_i");
                    outPortMap(padOp,"getNewData","getNewData_PadGen");
                    this->vhdl << Layer::getGetNewDataName(0) << " <= getNewData_PadGen;" << endl;
                }

                // finished
                this->outPortMap(padOp,"finished","finished_padding",false);
                DataGuard* datG = new DataGuard(target, 1, this->getSignalByName("featureCounter")->width(),(this->myArguments->getInputDepth()-1));
                this->addSubComponent(datG);
                this->inPortMap(datG,"Data_in","finished_padding");
                this->inPortMap(datG,"Guard_in","featureCounter");
                this->outPortMap(datG,"Data_out","finished",false);
                this->vhdl << instance(datG,"FinishedGuard");

                // newStep (also set newStepOrFinishedFeature = newStep or finished_padding)
                this->inPortMap(padOp,"newStep","newStepOrFinishedFeature");

                vhdl << declare("featureCounter_max",1) << " <= \"1\" when unsigned(featureCounter) = " << this->myArguments->getInputDepth()-1 << " else \"0\";" << endl;

                this->vhdl << "newStepOrFinishedFeature <= newStep or (finished_padding_rise and (not featureCounter_max));" << endl;

                // Prepare for Guards for validData_o-Signals
                this->outPortMap(padOp,"validData_o","validData_o_padding",true);

            }

            vhdl << instance(padOp,"paddingInstance_"+to_string(outFC));

            //pooling core
            PoolingCore* pooOp = new PoolingCore(target,myArguments->getWordSize(),myArguments->getCoreSize());

            addSubComponent(pooOp);
            for(int i=0; i<myArguments->getCoreSize()*myArguments->getCoreSize(); i++)
            {
                inPortMap(pooOp,"X"+to_string(i),"paddingOut_"+to_string(i)+"_"+to_string(outFC));
            }
            outPortMap(pooOp,"R","poolOut_"+to_string(outFC),true);
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


        if(myArguments->getActivationFunction()=="relu")
        {
            for(int outFC=0; outFC<forLimit; outFC++)
            {
                ReLU* relOp = new ReLU(target,this->myArguments->getWordSize());
                addSubComponent(relOp);
                inPortMap(relOp,"X","poolOut_"+to_string(outFC));
                outPortMap(relOp,"R","reluOut_"+to_string(outFC),true);
                vhdl << instance(relOp, "reluInstance_"+to_string(outFC));

                if(this->myArguments->getCalcAllParallel()==true)
                {
                    vhdl << "R_" << outFC << " <= reluOut_" << outFC << ";" << endl
                         << Layer::getValidDataName(outFC)+"_o" << " <= validData_temp;" << endl;
                }
                else
                {
                    if(this->outputMemoryParallelAccess==true)
                    {
                        // validData_o
                        for(unsigned int i=0;i<(unsigned int)this->myArguments->getInputDepth();i++)
                        {
                            DataGuard* datG = new DataGuard(target, 1, this->getSignalByName("featureCounter")->width(),i);
                            this->addSubComponent(datG);
                            this->inPortMap(datG,"Guard_in","featureCounter");
                            this->inPortMap(datG,"Data_in","validData_o_padding");
                            this->outPortMap(datG,"Data_out",Layer::getValidDataName(i)+"_o",false);
                            this->vhdl << instance(datG,"ValidData_o_DataGuard_"+to_string(i));
                        }
                        // Data out
                        for(unsigned int i=0;i<(unsigned int)this->myArguments->getInputDepth();i++)
                        {
                            vhdl << "R_" << i << " <= reluOut_" << outFC << ";" << endl;
                        }
                    }
                    else
                    {
                        this->vhdl << Layer::getValidDataName(0) << "_o <= validData_o_padding;" << endl;
                        this->vhdl << "R_0 <= reluOut_0;" << endl;
                    }
                }
            }
        }
        else
        {
            if(this->myArguments->getCalcAllParallel()==true)
            {
                for(int outFC=0; outFC<forLimit; outFC++)
                {
                    vhdl << "R_" << outFC << " <= poolOut_" << outFC << ";" << endl
                         << Layer::getValidDataName(outFC)+"_o" << " <= validData_temp;" << endl;
                }
            }
            else
            {
                if(this->outputMemoryParallelAccess==true)
                {
                    for(unsigned int i=0;i<(unsigned int)this->myArguments->getInputDepth();i++)
                    {
                        DataGuard* datG = new DataGuard(target, 1, this->getSignalByName("featureCounter")->width(),i);
                        this->addSubComponent(datG);
                        this->inPortMap(datG,"Guard_in","featureCounter");
                        this->inPortMap(datG,"Data_in","validData_o_padding");
                        this->outPortMap(datG,"Data_out",Layer::getValidDataName(i)+"_o",false);
                        this->vhdl << instance(datG,"ValidData_o_DataGuard_"+to_string(i));
                    }
                    for(unsigned int i=0;i<(unsigned int)this->myArguments->getInputDepth();i++)
                    {
                        vhdl << "R_" << i << " <= poolOut_0;" << endl;
                    }
                }
                else
                {
                    this->vhdl << Layer::getValidDataName(0) << "_o <= validData_o_padding;" << endl;
                    this->vhdl << "R_0 <= poolOut_0;" << endl;
                }
            }
        }

    }

    void PoolingLayer::createInputMemoryReset(Target *target)
    {
        int addressWidth = ceil(log2(this->myArguments->getInputHeight()*this->myArguments->getInputWidth()*this->myArguments->getInputDepth()));
        addOutput(Layer::getGetNewDataName(0)+"_reset",1);
        this->vhdl << Layer::getGetNewDataName(0) << "_reset <= counter_enable;" << endl;
        addOutput(Layer::getGetNewDataName(0)+"_reset_address",addressWidth);

        int counterWidth = 1;
        if(this->myArguments->getNumberOfOutputFeatures()>1) counterWidth = ceil(log2(this->myArguments->getNumberOfOutputFeatures()));

        // calculate the reset addresses
        map<unsigned int, unsigned int> addressMap;
        unsigned int featureSize = (unsigned int)(this->myArguments->getInputWidth()*this->myArguments->getInputHeight());
        for(int i=0; i<this->myArguments->getNumberOfOutputFeatures(); i++)
        {
            if(i<(this->myArguments->getNumberOfOutputFeatures()-1))
            {
                addressMap[i] = (i+1)*featureSize;
            }
            else
            {
                addressMap[i] = 0;
            }
        }

        // create lut
        GenericLut* addressLut = new GenericLut(target,"Address_Lut"+this->myArguments->getId(),addressMap,counterWidth,addressWidth);
        addSubComponent(addressLut);
        for(int i=0; i<counterWidth; i++)
        {
            this->vhdl << declare("counter_"+to_string(i)+"_std",1,false) << " <= featureCounter(" << i << ");" << endl;
            inPortMap(addressLut,"i"+to_string(i),"counter_"+to_string(i)+"_std");
        }
        for(int i=0; i<addressWidth; i++)
        {
            this->vhdl << Layer::getGetNewDataName(0) << "_reset_address(" << i << ") <= " << declare("reset_address_std_"+to_string(i)) << ";" << endl;
            outPortMap(addressLut,"o"+to_string(i),"reset_address_std_"+to_string(i), false);
        }
        this->vhdl << instance(addressLut,"Address_Lut_instance");
    }

}//namespace flopoco
