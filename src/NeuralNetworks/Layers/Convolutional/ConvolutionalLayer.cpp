// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of this Operator and the header of all nother ecessary Operators
#include "ConvolutionalLayer.hpp"

#include "NeuralNetworks/Layers/Convolutional/ConvolutionalCoreSimple.hpp"
#include "NeuralNetworks/Layers/Convolutional/ConvolutionalCoreWithWeightInputs.hpp"
#include "NeuralNetworks/Utility/WindowShiftRegister.hpp"
#include "NeuralNetworks/Utility/PaddingGenerator.hpp"
#include "NeuralNetworks/ActivationFunctions/ReLU.hpp"
#include "NeuralNetworks/Utility/Register.hpp"
#include "NeuralNetworks/Utility/BitheapWrapper.hpp"
#include "NeuralNetworks/Utility/ModuloCounter.hpp"
#include "NeuralNetworks/Utility/DataGuard.hpp"
#include "NeuralNetworks/MemoryManagement/DmaAccess/WeightFetcher.hpp"
#include "NeuralNetworks/Utility/Rounding.hpp"
#include "NeuralNetworks/NeuralNetwork.hpp"

#include "BitHeap/BitHeap.hpp"
#include "PrimitiveComponents/GenericMux.hpp"
#include "PrimitiveComponents/GenericLut.hpp"


using namespace std;
namespace flopoco {




    ConvolutionalLayer::ConvolutionalLayer(Target* target, NeuralNetwork* parent, int wordSize_, int fraction_, int horizontalSize_, int verticalSize_, int windowSize_, int numberOfOutputFeatures_, int numberOfInputFeatures_, vector <vector <vector <double> > > weights_, vector<double> biases_, int weightWordSize_, int weightFraction_, int paddingLeft_, string paddingType_, bool calcAllParallel_, string id_, int stride_, bool useAdderTree_, bool roundAfterConvCore_, char roundingType_, string activationFunction_, int paddingRight_, int paddingTop_, int paddingBot_, bool inputMemoryParallelAccess_, bool outputMemoryParallelAccess_) :
        Layer(target,parent), useAdderTree(useAdderTree_), roundAfterConvCore(roundAfterConvCore_), roundingType(roundingType_) {

        this->myArguments->setWordSize(wordSize_);
        this->myArguments->setFraction(fraction_);
        this->myArguments->setInputWidth(horizontalSize_);
        this->myArguments->setInputHeight(verticalSize_);
        this->myArguments->setNumberOfOutputFeatures(numberOfOutputFeatures_);
        this->myArguments->setInputDepth(numberOfInputFeatures_);
        this->myArguments->setWeightWordSize(weightWordSize_);
        this->myArguments->setWeightFraction(weightFraction_);
        this->myArguments->setCalcAllParallel(calcAllParallel_);
        this->myArguments->setActivationFunction(activationFunction_);
        this->myArguments->setLayerType("Convolutional");
        this->myArguments->setCoreSize(windowSize_);
        this->myArguments->setConvWeights(weights_);
        this->myArguments->setBiases(biases_);
        this->myArguments->setPaddingLeft(paddingLeft_);
        this->myArguments->setPaddingRight(paddingRight_);
        this->myArguments->setPaddingTop(paddingTop_);
        this->myArguments->setPaddingBot(paddingBot_);
        this->myArguments->setPaddingType(paddingType_);
        this->myArguments->setStride(stride_);
        this->myArguments->setId(id_);

        this->inputMemoryParallelAccess = inputMemoryParallelAccess_;
        this->outputMemoryParallelAccess = outputMemoryParallelAccess_;

        generateVHDLCode(target);
    }

    ConvolutionalLayer::ConvolutionalLayer(Target *target, NeuralNetwork* parent, LayerArguments *la, bool useAdderTree_, bool roundAfterConvCore_, char roundingType_, bool inputMemoryParallelAccess_, bool outputMemoryParallelAccess_) :
        Layer(target, la, parent), useAdderTree(useAdderTree_), roundAfterConvCore(roundAfterConvCore_), roundingType(roundingType_)
    {
        this->inputMemoryParallelAccess = inputMemoryParallelAccess_;
        this->outputMemoryParallelAccess = outputMemoryParallelAccess_;

        generateVHDLCode(target);
    }

    string ConvolutionalLayer::getOutputSignalName(int feature)
    {
        if(this->myArguments->getNumberOfOutputFeatures()<=feature)
        {
            cout << "ConvolutionalLayer.getOutputSignalName: requested port doesn't exist" << endl;
        }
        return "R"+to_string(feature);
    }

    string ConvolutionalLayer::getInputSignalName(int feature)
    {
        if(this->myArguments->getInputDepth()<=feature)
        {
            cout << "ConvolutionalLayer.getInputSignalName: requested port doesn't exist" << endl;
        }
        return "X"+to_string(feature);
    }

    string ConvolutionalLayer::getIntermediateResultName(unsigned int featureNumber)
    {
        if((unsigned int)this->myArguments->getInputDepth()<=featureNumber)
        {
            cout << "ConvolutionalLayer.getIntermediateResultName: requested port doesn't exist" << endl;
        }
        return this->getOutputSignalName(featureNumber)+"_intermediate";
    }

    void ConvolutionalLayer::buildShiftAndPadParallel(Target* target)
    {

        for(int featureCounter=0; featureCounter<myArguments->getInputDepth(); featureCounter++)
        {
            // declare inputs
            addInput("X"+to_string(featureCounter), myArguments->getWordSize());
            addInput(Layer::getValidDataName(featureCounter)+"_i",1);
            addOutput(Layer::getGetNewDataName(featureCounter),1);

            // shift register
            WindowShiftRegister* shiftReg = new WindowShiftRegister(target,myArguments->getWordSize(),myArguments->getCoreSize(),myArguments->getInputWidth());
            addSubComponent(shiftReg);
            inPortMap(shiftReg,"X","X"+to_string(featureCounter));
            inPortMap(shiftReg,"enable",Layer::getValidDataName(featureCounter)+"_i");
            inPortMap(shiftReg,"newStep","newStep");
            for(int portC=0;portC<myArguments->getCoreSize()*myArguments->getCoreSize();portC++)
            {
                outPortMap(shiftReg,"R"+to_string(portC),"PaddingGenInput_inputFeature_"+to_string(featureCounter)+"_number_"+to_string(portC),true);
            }
            vhdl << instance(shiftReg,"WindowShiftRegister"+to_string(featureCounter)+"_instance") << endl;

            // padding generator
            PaddingGenerator* padGen = new PaddingGenerator(target,myArguments->getWordSize(),myArguments->getCoreSize(),myArguments->getInputWidth(),myArguments->getInputHeight(),myArguments->getPaddingTop(),myArguments->getStride(),myArguments->getPaddingType(),(featureCounter==0),false,myArguments->getPaddingBot(),myArguments->getPaddingLeft(),myArguments->getPaddingRight(),(featureCounter==0?true:false));
            addSubComponent(padGen);
            for(int portC=0; portC<myArguments->getCoreSize()*myArguments->getCoreSize(); portC++)
            {
                inPortMap(padGen,"X"+to_string(portC),"PaddingGenInput_inputFeature_"+to_string(featureCounter)+"_number_"+to_string(portC));
                outPortMap(padGen,"R"+to_string(portC),"PaddingGenOutput_inputFeature_"+to_string(featureCounter)+"_number_"+to_string(portC),true);
            }
            inPortMap(padGen,"validData_i",Layer::getValidDataName(featureCounter)+"_i");
            inPortMap(padGen,"newStep","newStep");
            outPortMap(padGen,"getNewData",Layer::getGetNewDataName(featureCounter),false);
            if(featureCounter==0)
            {
                outPortMap(padGen,"validData_o","validData_temp",true);
                outPortMap(padGen,"finished","finished_tmp",true);
            }
            vhdl << instance(padGen,"PaddingGen_"+to_string(featureCounter)+"_instance");

        }
    }

    void ConvolutionalLayer::buildConvCoresParallel(Target* target)
    {
        int convCoreIdCounter=0;
        for(int featureCounter=0; featureCounter<myArguments->getInputDepth(); featureCounter++)
        {
            // one ConvCore for each OUTPUT feature
            // (get pipelinestages from all ConvCores for each OUTPUT-feature and insert additional Delays if needed)
            vector <ConvolutionalCoreSimple*> tempConvCores;
            this->ConvCores.push_back(tempConvCores);
            int tempOutWordSize;
            int tempOutFraction;
            for(int outFeatureCounter=0; outFeatureCounter<myArguments->getNumberOfOutputFeatures(); outFeatureCounter++)
            {
                ConvolutionalCoreSimple* convCore = new ConvolutionalCoreSimple(target,myArguments->getWordSize(),myArguments->getFraction(),myArguments->getWeightWordSize(),myArguments->getWeightFraction(),myArguments->getCoreSize(),myArguments->getConvWeights()[featureCounter][outFeatureCounter],useAdderTree,(roundAfterConvCore==true?roundingType:0xFF),this->myArguments->getId()+"_"+to_string(convCoreIdCounter));
                convCoreIdCounter++;
                this->addSubComponent(convCore);
                for(int portC=0; portC<myArguments->getCoreSize()*myArguments->getCoreSize(); portC++)
                {
                    this->inPortMap(convCore,"X"+to_string(portC),"PaddingGenOutput_inputFeature_"+to_string(featureCounter)+"_number_"+to_string(portC));
                }
                this->outPortMap(convCore,"R","ConvCoreOutput_inputFeature_"+to_string(featureCounter)+"_outputFeature_"+to_string(outFeatureCounter),true);
                this->vhdl << this->instance(convCore,"ConvCore_inputFeature_"+to_string(featureCounter)+"_outputFeature_"+to_string(outFeatureCounter)+"_instance");
                this->ConvCores[featureCounter].push_back(convCore);

                if(outFeatureCounter>0)
                {
                    if(tempOutWordSize!=(int)convCore->getOutputWordSize() || tempOutFraction!=(int)convCore->getOutputFraction())
                    {
                        stringstream e;
                        e << "The Convolutional Cores have a different outputWordSize!";
                        THROWERROR(e.str());
                    }
                }
                tempOutWordSize=convCore->getOutputWordSize();
                tempOutFraction=convCore->getOutputFraction();
            }
        }

        // delay data and flag-signals for all convCores with pipeline stages < maxPipeline stages
        unsigned int inputFCCounter=0;
        unsigned int outputFCCounter=0;
        for(unsigned int inputFC=0; inputFC<this->ConvCores.size(); inputFC++)
        {
            for(unsigned int outputFC=0; outputFC<this->ConvCores[inputFC].size(); outputFC++)
            {
                if(ConvCores[inputFC][outputFC]->getPipelineDepth()>=ConvCores[inputFCCounter][outputFCCounter]->getPipelineDepth())
                {
                    inputFCCounter=inputFC;
                    outputFCCounter=outputFC;
                }
            }
        }
        string tempSignalName="ConvCoreOutput_inputFeature_"+to_string(inputFCCounter)+"_outputFeature_"+to_string(outputFCCounter); // this is the output-signal from the ConvCore with the most internal pipelinestages
        this->syncCycleFromSignal(tempSignalName);
    }

    void ConvolutionalLayer::buildAdderTreesParallel(Target* target)
    {
        for(int outFeatureCounter=0; outFeatureCounter<this->myArguments->getNumberOfOutputFeatures(); outFeatureCounter++)
        {
            addOutput(Layer::getValidDataName(outFeatureCounter)+"_o",1);

            // with bitheap wrapper
            vector<int> wIns;
            vector<bool> signs;
            vector<int> weights;
            for(int i=0; i<this->myArguments->getInputDepth(); i++)
            {
                int wIn = this->ConvCores[i][outFeatureCounter]->getOutputWordSize();
                wIns.push_back(wIn);
                signs.push_back(true);
                weights.push_back(0);
            }
            wIns.push_back(this->myArguments->getWeightWordSize());
            signs.push_back(true);
            weights.push_back(this->myArguments->getFraction());
            BitheapWrapper* bitH = new BitheapWrapper(target,wIns,this->myArguments->getInputDepth()+1,signs,weights);
            addSubComponent(bitH);
            // connect conv cores
            for(int i=0; i<this->myArguments->getInputDepth(); i++)
            {
                inPortMap(bitH,"X"+to_string(i),"ConvCoreOutput_inputFeature_"+to_string(i)+"_outputFeature_"+to_string(outFeatureCounter));
            }
            // create and connect bias
            double scal = this->myArguments->getBias(outFeatureCounter) * pow(2,this->myArguments->getWeightFraction());
            this->vhdl << declare("Bias_"+to_string(outFeatureCounter),this->myArguments->getWeightWordSize()) << " <= std_logic_vector(to_signed(" << (int)round(scal) << "," << this->myArguments->getWeightWordSize() << "));" << endl;
            inPortMap(bitH,"X"+to_string(this->myArguments->getInputDepth()),"Bias_"+to_string(outFeatureCounter));
            // connect output
            outPortMap(bitH,"R","OutputFeature_"+to_string(outFeatureCounter)+"_Result",true);
            // create instance in vhdl
            this->vhdl << instance(bitH,"Bitheap_"+to_string(outFeatureCounter));
         }
        //sync
        syncCycleFromSignal("OutputFeature_0_Result");
        //round
        for(int outFeatureCounter=0; outFeatureCounter<this->myArguments->getNumberOfOutputFeatures(); outFeatureCounter++)
        {
            roundOutput(target,this->getSignalByName("OutputFeature_"+to_string(outFeatureCounter)+"_Result")->width(),this->ConvCores[0][outFeatureCounter]->getOutputFraction(),"OutputFeature_"+to_string(outFeatureCounter)+"_Result","feature"+to_string(outFeatureCounter)+"_rounded","Rounding_instance_"+to_string(outFeatureCounter),this->roundingType);
        }
    }

    void ConvolutionalLayer::roundOutput(Target* target, int wordSizeFrom, int fractionFrom, string signalNameFrom, string signalNameTo, string instanceName, char round)
    {
        flopoco::roundingTypeEnum rType;
        if(roundingType==0x00)
        {
            rType = flopoco::roundingTypeEnum::truncation;
        }
        else if(roundingType==0x01)
        {
            rType = flopoco::roundingTypeEnum::saturation;
        }
        else
        {
            cout << "ConvolutionalCoreSimple.roundOutput: atm only truncation/saturation are implemented. FloPoCo is rounding with truncation instead" << endl;
            rType = flopoco::roundingTypeEnum::truncation;
        }
        Rounding* roundOp = new Rounding(target,wordSizeFrom,fractionFrom,this->myArguments->getWordSize(),this->myArguments->getFraction(),rType);
        addSubComponent(roundOp);
        inPortMap(roundOp,"X",signalNameFrom);
        outPortMap(roundOp,"R",signalNameTo,true);
        this->vhdl << instance(roundOp,instanceName);

    }

    void ConvolutionalLayer::buildActivationFunction(Target* target, int outFeatureCounter)
    {
        if(this->myArguments->getActivationFunction()=="relu")
        {
            ReLU* relu = new ReLU(target, this->myArguments->getWordSize());
            addSubComponent(relu);
            inPortMap(relu,"X","feature"+to_string(outFeatureCounter)+"_rounded");
            outPortMap(relu,"R","ReLU"+to_string(outFeatureCounter)+"_out",true);
            vhdl << instance(relu,"ReLU"+to_string(outFeatureCounter)+"_instance");

            addOutput("R"+to_string(outFeatureCounter),this->myArguments->getWordSize());
            vhdl << "R" << outFeatureCounter << " <= ReLU" << outFeatureCounter << "_out;" << endl;
        }
        else
        {
            addOutput("R"+to_string(outFeatureCounter),this->myArguments->getWordSize());
            vhdl << "R" << outFeatureCounter << " <= feature" << outFeatureCounter << "_rounded;" << endl;
        }
    }

    void ConvolutionalLayer::buildParallelCalculation(Target* target)
    {
        if(this->inputMemoryParallelAccess==false || this->outputMemoryParallelAccess==false)
        {
            THROWERROR("Cant build Parallel Calculation, if input/output features can't be accessed parallel");
        }
        // build one convolutional chain for each INPUT feature!
        buildShiftAndPadParallel(target);
        buildConvCoresParallel(target);

        // Add all ConvCore-outputs for each OUTPUT-feature up and round back to original word size (while doing that, also create validData_o-signals because we need them anyways)
        this->buildAdderTreesParallel(target);

        // delay finished-signal and validData-signal for x cycles, with x=maxNumberOfPipelineStages(all AdderTrees) and assign outputs
        vhdl << "finished <= finished_tmp;" << endl;

        for(int featureCounter=0; featureCounter<myArguments->getNumberOfOutputFeatures(); featureCounter++)
        {
            vhdl << Layer::getValidDataName(featureCounter)+"_o" << " <= validData_temp;" << endl;
            this->buildActivationFunction(target, featureCounter);
        }
    }

    void ConvolutionalLayer::buildSerialCalculation(Target *target)
    {
        // ports
        this->createPortsSerial(target);
        // signals
        this->declareSignals(target);
        // data guards and multiplexers for BRAM communication
        this->createDataGuardsAndMultiplexersSerial(target);
        // convolutional core
        ConvolutionalCoreWithWeightInputs* conv = this->createConvCoreSerial(target);
        // biases
        this->createBiasesSerial(target);
        // bit heap
        unsigned int pipelineStagesBitH = this->createBitHeapSerial(target, conv->getOutputWordSize());
        // pipeline registers
        this->createPipelineRegistersSerial(target,conv->getPipelineDepth(),pipelineStagesBitH);
        // shift and pad
        this->createShitAndPadSerial(target);
        // DMA controller
        this->createDMAControllerSerial(target);
        // feature counters
        this->createCountersSerial(target);
        // activation function
        this->createActivationFunctionSerial(target);
        // LUT for intermediate data resetting
        this->createLutForIntermediateReset(target);
    }

    void ConvolutionalLayer::createPortsSerial(Target* target)
    {
        // newStep and finished already exist... => finished must still be assigned
        this->vhdl << "finished <= finished_delayed and " << ConvolutionalLayer::getOutputFeatureCounterOutputSignalName() << "_Equals_Max and " << ConvolutionalLayer::getInputFeatureCounterOutputSignalName() << "_Equals_Max;" << endl;

        if(this->inputMemoryParallelAccess==false)
        {
            // There is ONE memory controller, that stores data for EVERY input feature
            // X
            addInput(this->getInputSignalName(0),this->myArguments->getWordSize());
            // validData_i
            addInput(Layer::getValidDataName(0)+"_i",1,true);
            // getNewData
            addOutput(Layer::getGetNewDataName(0),1,true);
            // getNewData_reset
            addOutput(Layer::getGetNewDataName(0)+"_reset",1,true);
            // getNewData_reset_address => this is needed because the padding every time reads more than it should read and ignores it (because it's easier that way), so after each feature, the correct input address has to be set
            addOutput(Layer::getGetNewDataName(0)+"_reset_address",(ceil(log2(this->myArguments->getInputWidth()*this->myArguments->getInputHeight()*this->myArguments->getInputDepth()))),true);
        }
        else
        {
            // each input feature has its own memory controller
            for(int i=0;i<this->myArguments->getInputDepth();i++)
            {
                // X
                addInput(this->getInputSignalName(i),this->myArguments->getWordSize());
                // validData_i
                addInput(Layer::getValidDataName(i)+"_i",1,true);
                // getNewData
                addOutput(Layer::getGetNewDataName(i),1,true);
                // getNewData_reset
                addOutput(Layer::getGetNewDataName(i)+"_reset",1,true);
            }
        }


        if(this->outputMemoryParallelAccess==false)
        {
            // There is ONE memory controller, that stores data for EVERY output feature
            // R
            addOutput(this->getOutputSignalName(0),this->myArguments->getWordSize());
            // validData_o
            addOutput(Layer::getValidDataName(0)+"_o",1,true);
            // getNewDataIntermediate
            addOutput(Layer::getGetNewDataName(0)+"_intermediate",1,true);
            // getNewDataIntermediate_reset
            addOutput(Layer::getGetNewDataName(0)+"_intermediate_reset",1,true);
            // R_intermediate
            addInput(this->getOutputSignalName(0)+"_intermediate",this->myArguments->getWordSize(),true);
            // reset_address
            this->outputWidth = NeuralNetwork::calculateNewOutputSize(this->myArguments->getInputWidth(),this->myArguments->getPaddingLeft(), this->myArguments->getPaddingRight(), this->myArguments->getCoreSize(), this->myArguments->getStride());
            this->outputHeight = NeuralNetwork::calculateNewOutputSize(this->myArguments->getInputHeight(),this->myArguments->getPaddingTop(), this->myArguments->getPaddingBot(), this->myArguments->getCoreSize(), this->myArguments->getStride());
            this->addressWidth = ceil(log2(this->outputWidth * this->outputHeight * this->myArguments->getNumberOfOutputFeatures()));
            addOutput("reset_address",this->addressWidth);
        }
        else
        {
            this->addressWidth = 0; // just set addressWidth to 0, because no reset_address is needed
            // each output feature has its own memory controller
            for(int i=0;i<this->myArguments->getNumberOfOutputFeatures();i++)
            {
                // R
                addOutput(this->getOutputSignalName(i),this->myArguments->getWordSize());
                // validData_o
                addOutput(Layer::getValidDataName(i)+"_o",1,true);
                // getNewDataIntermediate
                addOutput(Layer::getGetNewDataName(i)+"_intermediate",1,true);
                // getNewDataIntermediate_reset
                addOutput(Layer::getGetNewDataName(i)+"_intermediate_reset",1,true);
                // R_intermediate
                addInput(this->getOutputSignalName(i)+"_intermediate",this->myArguments->getWordSize(),true);
            }
        }

        // to connect Weight fetcher with global DMA controller
        // weight input
        addInput(Layer::getWeightInputName(),32);
        // valid data input
        addInput(Layer::getWeightsValidInputName(),1);
        // last weights input
        addInput(Layer::getLastWeightsInputName(),1);
        // start new read access output
        addOutput(Layer::getNewReadAccessOutputName(),1);
        // number of bytes output
        addOutput(Layer::getNumberOfBytesOutputName(),23);
        // start address output
        addOutput(Layer::getNextStartAddressOutputName(),32);
    }

    void ConvolutionalLayer::declareSignals(Target* target)
    {
        this->declare(ConvolutionalLayer::getOutputFeatureCounterOutputSignalName(),this->getOutputCounterWidth());
        this->declare(ConvolutionalLayer::getInputFeatureCounterOutputSignalName(),this->getInputCounterWidth());
        declare("getNewData",1);
        declare("X",this->myArguments->getWordSize());
        declare("validData_i",1);
        declare("Intermediate_Result",this->myArguments->getWordSize());
        declare("Conv_Core_Weight_Enable",1);
        declare("Bias",this->myArguments->getWeightWordSize());
        declare("Second_Input_Mux_Select",1);
        declare("finished_delayed",1);
        declare("validData_delayed_conv",1);
        declare("validData_delayed",1);
        declare("Read_New_Weights_temp",1);
        declare("finished_temp",1);
        declare("validData_o_tmp",1);
        declare("Read_New_Weights",1);
        declare("newStep_shiftAndPad",1);
        declare("New_Weights_Ready",1);
        declare(ConvolutionalLayer::getOutputFeatureCounterOutputSignalName()+"_Equals_Max",1);
        declare(ConvolutionalLayer::getInputFeatureCounterOutputSignalName()+"_Equals_Max",1);
        declare("finished_delayed_fall",1);
        declare("counter_reset",1);
        declare("Enable_Output_Feature_Counter",1);
        declare("R_temp",this->myArguments->getWordSize());
    }

    void ConvolutionalLayer::createDataGuardsAndMultiplexersSerial(Target* target)
    {
        ////////////////
        // input side //
        ////////////////
        if(this->inputMemoryParallelAccess==false)
        {
            this->vhdl << "X <= " << this->getInputSignalName(0) << ";" << endl;
            this->vhdl << "validData_i <= " << Layer::getValidDataName(0) << "_i;" << endl;
            this->vhdl << Layer::getGetNewDataName(0) << " <= getNewData;" << endl;
        }
        else
        {
            if(this->myArguments->getInputDepth() > 1)
            {
                // Multiplexers to select the right input/validData_i
                GenericMux* inputMux = new GenericMux(target,this->myArguments->getWordSize(),this->myArguments->getInputDepth());
                GenericMux* validIMux = new GenericMux(target,1,this->myArguments->getInputDepth());
                addSubComponent(inputMux);
                addSubComponent(validIMux);
                // data inputs
                for(int i=0;i<this->myArguments->getInputDepth();i++)
                {
                    inPortMap(inputMux,inputMux->getInputName(i),this->getInputSignalName(i));
                    inPortMap(validIMux,validIMux->getInputName(i), Layer::getValidDataName(i)+"_i");
                }
                // select input
                inPortMap(inputMux,inputMux->getSelectName(),ConvolutionalLayer::getInputFeatureCounterOutputSignalName());
                inPortMap(validIMux,validIMux->getSelectName(),ConvolutionalLayer::getInputFeatureCounterOutputSignalName());
                // data output
                outPortMap(inputMux,inputMux->getOutputName(),"X",false);
                outPortMap(validIMux,validIMux->getOutputName(),"validData_i",false);

                this->vhdl << instance(inputMux,"Input_Multiplexer");
                this->vhdl << instance(validIMux,"Valid_Data_Multiplexer");
            }
            else
            {
                this->vhdl << "X <= " << this->getInputSignalName(0) << ";" << endl;
                this->vhdl << "validData_i <= " << Layer::getValidDataName(0) << "_i;" << endl;
            }

            // Data Guards for getNewData (from input) signals
            for(int i=0; i<this->myArguments->getInputDepth(); i++)
            {
                DataGuard* dat = new DataGuard(target,1,this->getInputCounterWidth(),i);
                addSubComponent(dat);

                inPortMap(dat,"Data_in","getNewData");
                inPortMap(dat,"Guard_in",ConvolutionalLayer::getInputFeatureCounterOutputSignalName());
                outPortMap(dat,"Data_out",Layer::getGetNewDataName(i),false);
                this->vhdl << instance(dat,"Data_Guard_GetNewData_"+to_string(i));
            }
        }

        /////////////////
        // output side //
        /////////////////
        if(this->outputMemoryParallelAccess==false)
        {
            this->vhdl << Layer::getValidDataName(0) << "_o <= validData_delayed;" << endl;
            this->vhdl <<  Layer::getGetNewDataName(0) << "_intermediate <= validData_delayed_conv;" << endl;
            this->vhdl << "Intermediate_Result <= " << this->getOutputSignalName(0)+"_intermediate;" << endl;
        }
        else
        {
            // Data Guards for validData_o/getNewData_intermediate (for outputs) signals
            for(int i=0; i<this->myArguments->getNumberOfOutputFeatures(); i++)
            {
                // validData_o
                DataGuard* dat = new DataGuard(target,1,this->getOutputCounterWidth(),i);
                addSubComponent(dat);
                inPortMap(dat,"Data_in","validData_delayed");
                inPortMap(dat,"Guard_in",ConvolutionalLayer::getOutputFeatureCounterOutputSignalName());
                outPortMap(dat,"Data_out",Layer::getValidDataName(i)+"_o",false);
                this->vhdl << instance(dat,"Data_Guard_ValidDataO_"+to_string(i));

                // getNewData_intermediate
                DataGuard* dat2 = new DataGuard(target,1,this->getOutputCounterWidth(),i);
                addSubComponent(dat2);
                inPortMap(dat2,"Data_in","validData_delayed_conv");
                inPortMap(dat2,"Guard_in",ConvolutionalLayer::getOutputFeatureCounterOutputSignalName());
                outPortMap(dat2,"Data_out",Layer::getGetNewDataName(i)+"_intermediate",false);
                this->vhdl << instance(dat2,"Data_Guard_"+Layer::getGetNewDataName(i)+"_intermediate_"+to_string(i));
            }

            if(this->myArguments->getNumberOfOutputFeatures() > 1)
            {
                // mux for intermediate result
                GenericMux* intermediateMux = new GenericMux(target,this->myArguments->getWordSize(),this->myArguments->getNumberOfOutputFeatures());
                addSubComponent(intermediateMux);
                inPortMap(intermediateMux,intermediateMux->getSelectName(),ConvolutionalLayer::getOutputFeatureCounterOutputSignalName());
                outPortMap(intermediateMux,intermediateMux->getOutputName(),"Intermediate_Result",false);
                for(int i=0; i<this->myArguments->getNumberOfOutputFeatures(); i++)
                {
                    inPortMap(intermediateMux,intermediateMux->getInputName(i),this->getOutputSignalName(i)+"_intermediate");
                }
                this->vhdl << instance(intermediateMux,"Intermediate_Result_Mux");
            }
            else
            {
                this->vhdl << "Intermediate_Result <= " << this->getOutputSignalName(0)+"_intermediate;" << endl;
            }
        }

    }

    ConvolutionalCoreWithWeightInputs* ConvolutionalLayer::createConvCoreSerial(Target* target)
    {
        ConvolutionalCoreWithWeightInputs* conv = new ConvolutionalCoreWithWeightInputs(target,this->parent->getMultiplicationType(),this->myArguments->getWordSize(),this->myArguments->getFraction(),this->myArguments->getWeightWordSize(),this->myArguments->getWeightFraction(),this->myArguments->getCoreSize(),((this->roundAfterConvCore==true)?(this->roundingType):(0xFF)),this->myArguments->getId());
        addSubComponent(conv);
        for(int i = 0; i < this->myArguments->getCoreSize()*this->myArguments->getCoreSize(); i++)
        {
            string portName = getConvCoreInputSignalName(i);
            declare(portName,this->myArguments->getWordSize());
            inPortMap(conv,"X"+to_string(i),portName);

            portName = ConvolutionalLayer::getConvCoreWeightSignalName(i);
            declare(portName,this->myArguments->getWeightWordSize());
            inPortMap(conv,"W"+to_string(i),portName);
        }
        declare(ConvolutionalLayer::getConvCoreOutputSignalName(),conv->getOutputWordSize());
        outPortMap(conv,"R",ConvolutionalLayer::getConvCoreOutputSignalName(),false);
        inPortMap(conv,"Weight_Enable","Conv_Core_Weight_Enable");
        this->vhdl << instance(conv,"ConvolutionalCore");

        return conv;
    }

    void ConvolutionalLayer::createBiasesSerial(Target* target)
    {
        // copy of bias vector
        vector<double> biases = this->myArguments->getBiases();

        if(biases.size()==0) THROWERROR("Bias Vector can't be empty!");
        if(biases.size()>1)
        {
            // create lut to select the right one
            map<unsigned int, unsigned int> biasMap;
            for(unsigned int i=0; i<biases.size(); i++)
            {
                double scaledBias = biases[i] * pow(2,this->myArguments->getWeightFraction());
                if((double)abs(scaledBias) >= pow(2,32))
                {
                    THROWERROR("Bias " + to_string(i) + " doesn't fit word size");
                }
                biasMap[i] = (unsigned int)((int)scaledBias); // looks stupid, but works
            }
            GenericLut* biasLut = new GenericLut(target,"Bias_Lut"+this->myArguments->getId(),biasMap,this->getOutputCounterWidth(),this->myArguments->getWeightWordSize());
            addSubComponent(biasLut);
            // in port map
            for (unsigned int i = 0; i < this->getOutputCounterWidth(); i++) {
                this->vhdl << this->declare("Bias_Lut_i" + to_string(i), 1, false) << " <= "
                           << ConvolutionalLayer::getOutputFeatureCounterOutputSignalName() << "(" << i << ");" << endl;
                inPortMap(biasLut, "i" + to_string(i), "Bias_Lut_i" + to_string(i));
            }
            // out port map
            for (unsigned int i = 0; i < this->myArguments->getWeightWordSize(); i++)
            {
                outPortMap(biasLut,"o"+to_string(i),"Bias_"+to_string(i)+"_std",true);
                this->vhdl << "Bias(" << i << ") <= Bias_" << i << "_std;" << endl;
            }
            this->vhdl << instance(biasLut,"Bias_Lut_instance");
        }
        else
        {
            double scaledBias = biases[0] * pow(2,this->myArguments->getWeightFraction());
            this->vhdl << "Bias <= std_logic_vector(to_signed(" << ((int)scaledBias) << "," << this->myArguments->getWeightWordSize() << "));" << endl;
        }
    }

    unsigned int ConvolutionalLayer::createBitHeapSerial(Target* target, unsigned int wInConv)
    {
        // bring bias and intermediate result on the same word size
        // get MSB and LSB
        int fractionBias = this->myArguments->getWeightFraction();
        int fractionInte = this->myArguments->getFraction();
        int biggerFraction = fractionInte;
        if(fractionBias > biggerFraction) biggerFraction = fractionBias;
        int msbBias = this->myArguments->getWeightWordSize() - fractionBias;
        int msbInte = this->myArguments->getWordSize() - fractionInte;
        int biggerMSB = msbInte;
        if(msbBias > biggerMSB) biggerMSB = msbBias;
        int newWordSize = biggerFraction + biggerMSB;

        // resize bias
        this->vhdl << declare("Bias_resized",newWordSize) << " <= ";
        if(msbBias < biggerMSB)
        {
            for(int i=0; i<biggerMSB-msbBias; i++)
            {
                this->vhdl << "Bias(" << this->myArguments->getWeightWordSize()-1 << ") & ";
            }
        }
        this->vhdl << "Bias";
        if(fractionBias < biggerFraction)
        {
            this->vhdl << " & \"";
            for(int i=0; i<biggerFraction-fractionBias; i++)
            {
                this->vhdl << "0";
            }
            this->vhdl << "\"";
        }
        this->vhdl << ";" << endl;

        // resize intermediate result
        this->vhdl << declare("Intermediate_Result_resized",newWordSize) << " <= ";
        if(msbInte < biggerMSB)
        {
            for(int i=0; i<biggerMSB-msbInte; i++)
            {
                this->vhdl << "Intermediate_Result(" << this->myArguments->getWordSize()-1 << ") & ";
            }
        }
        this->vhdl << "Intermediate_Result";
        if(fractionInte < biggerFraction)
        {
            this->vhdl << " & \"";
            for(int i=0; i<biggerFraction-fractionInte; i++)
            {
                this->vhdl << "0";
            }
            this->vhdl << "\"";
        }
        this->vhdl << ";" << endl;

        // check if the input feature counter equals zero (select signal for the upcoming mux)
        this->vhdl << "process(" << ConvolutionalLayer::getInputFeatureCounterOutputSignalName() << ")" << endl;
        this->vhdl << "begin" << endl;
        this->vhdl << tab << "if(unsigned(" << ConvolutionalLayer::getInputFeatureCounterOutputSignalName() << ") = 0) then" << endl;
        this->vhdl << tab << tab << "Second_Input_Mux_Select <= \"1\";" << endl;
        this->vhdl << tab << "else" << endl;
        this->vhdl << tab << tab << "Second_Input_Mux_Select <= \"0\";" << endl;
        this->vhdl << tab << "end if;" << endl;
        this->vhdl << "end process;" << endl;

        // create mux to select second input of bitheap
        GenericMux* secondInputMux = new GenericMux(target,newWordSize,2);
        addSubComponent(secondInputMux);
        inPortMap(secondInputMux,secondInputMux->getSelectName(),"Second_Input_Mux_Select");
        inPortMap(secondInputMux,secondInputMux->getInputName(1),"Bias_resized");
        inPortMap(secondInputMux,secondInputMux->getInputName(0),"Intermediate_Result_resized");
        outPortMap(secondInputMux,secondInputMux->getOutputName(),"Second_Input_Bitheap",true);
        this->vhdl << instance(secondInputMux,"Second_Input_Mux");

        // create bit heap
        vector<bool> signs;
        signs.push_back(true);
        signs.push_back(true);
        vector<int> wIn;
        wIn.push_back(wInConv);
        wIn.push_back(newWordSize);
        vector<int> weights;
        weights.push_back(0);
        weights.push_back(fractionBias+fractionInte-biggerFraction);
        BitheapWrapper* bitH = new BitheapWrapper(target,wIn,2,signs,weights);
        addSubComponent(bitH);
        inPortMap(bitH,"X0",ConvolutionalLayer::getConvCoreOutputSignalName());
        inPortMap(bitH,"X1","Second_Input_Bitheap");
        declare("Bitheap_out",bitH->getOutputWordSize());
        outPortMap(bitH,"R","Bitheap_out",false);
        this->vhdl << instance(bitH,"Bitheap_Wrapper");
        return bitH->getPipelineDepth();
    }

    void ConvolutionalLayer::createPipelineRegistersSerial(Target* target, unsigned int pipelineStagesConv, unsigned int pipelineStagesBitH)
    {
        if(pipelineStagesConv<1)
        {
            THROWERROR("pipeline stages of conv core must be >= 1 for this convolutional layer architecture to work (try increasing the target frequency)");
        }
        // validData_delayed_conv already exists
        // validData_delayed already exists

        Register* readNewWeightsReg = new Register(target,1,pipelineStagesBitH+pipelineStagesConv+1,false);
        Register* finishedReg = new Register(target,1,pipelineStagesBitH+pipelineStagesConv);
        Register* validDataRegConv = new Register(target,1,pipelineStagesConv-1);
        Register* validDataRegTotal = new Register(target,1,pipelineStagesBitH+1);

        addSubComponent(readNewWeightsReg);
        addSubComponent(finishedReg);
        addSubComponent(validDataRegConv);
        addSubComponent(validDataRegTotal);

        inPortMap(readNewWeightsReg,"X","Read_New_Weights_temp");
        inPortMap(finishedReg,"X","finished_temp");
        inPortMap(validDataRegConv,"X","validData_o_tmp");
        inPortMap(validDataRegTotal,"X","validData_delayed_conv");

        outPortMap(readNewWeightsReg,"R","Read_New_Weights",false);
        outPortMap(finishedReg,"R","finished_delayed",false);
        outPortMap(validDataRegConv,"R","validData_delayed_conv",false);
        outPortMap(validDataRegTotal,"R","validData_delayed",false);

        this->vhdl << instance(readNewWeightsReg,"Read_New_Weights_Reg");
        this->vhdl << instance(finishedReg,"Finished_Reg");
        this->vhdl << instance(validDataRegConv,"ValidData_Reg_Conv");
        this->vhdl << instance(validDataRegTotal,"ValidData_Reg_Total");
    }

    void ConvolutionalLayer::createShitAndPadSerial(Target* target)
    {
        // reset intermediate result address
        if(this->outputMemoryParallelAccess==true)
        {
            for(int i=0; i<this->myArguments->getNumberOfOutputFeatures(); i++)
            {
                this->vhdl << Layer::getGetNewDataName(i) << "_intermediate_reset <= newStep_shiftAndPad;" << endl;
            }
        }
        else
        {
            this->vhdl << Layer::getGetNewDataName(0) << "_intermediate_reset <= newStep_shiftAndPad and (not " << ConvolutionalLayer::getInputFeatureCounterOutputSignalName() << "_Equals_Max);" << endl;
        }

        // Shift Register
        unsigned int numberOfInputs = this->myArguments->getCoreSize()*this->myArguments->getCoreSize();
        WindowShiftRegister* shift = new WindowShiftRegister(target,this->myArguments->getWordSize(),this->myArguments->getCoreSize(),this->myArguments->getInputWidth());
        addSubComponent(shift);
        inPortMap(shift,"X","X");
        inPortMap(shift,"enable","validData_i");
        inPortMap(shift,"newStep","newStep_shiftAndPad");
        for(unsigned int i=0; i<numberOfInputs; i++)
        {
            outPortMap(shift,"R"+to_string(i),declare("R_shift_"+to_string(i),this->myArguments->getWordSize()),false);
        }
        this->vhdl << instance(shift,"Shift_Register");

        // Padding
        PaddingGenerator* pad = new PaddingGenerator(target,this->myArguments->getWordSize(),this->myArguments->getCoreSize(),
                this->myArguments->getInputWidth(),this->myArguments->getInputHeight(),this->myArguments->getPaddingTop(),
                this->myArguments->getStride(),this->myArguments->getPaddingType(),true,true,this->myArguments->getPaddingBot(),
                this->myArguments->getPaddingLeft(),this->myArguments->getPaddingRight(),this->myArguments->getStride());
        addSubComponent(pad);
        for(unsigned int i=0; i<numberOfInputs; i++)
        {
            inPortMap(pad,"X"+to_string(i),"R_shift_"+to_string(i));
        }
        inPortMap(pad,"newWeightsReady","New_Weights_Ready");
        inPortMap(pad,"validData_i","validData_i");
        inPortMap(pad,"newStep","newStep_shiftAndPad");

        for(unsigned int i=0; i<numberOfInputs; i++)
        {
            outPortMap(pad,"R"+to_string(i),ConvolutionalLayer::getConvCoreInputSignalName(i),false);
        }
        outPortMap(pad,"readNewWeights","Read_New_Weights_temp",false);
        outPortMap(pad,"convCoreWeightEnable","Conv_Core_Weight_Enable",false);
        outPortMap(pad,"getNewData","getNewData",false);
        outPortMap(pad,"finished","finished_temp",false);
        outPortMap(pad,"validData_o","validData_o_tmp",false);
        this->vhdl << instance(pad,"Padding");

        // assign newStep signal
        this->vhdl << "newStep_shiftAndPad <= newStep or (finished_delayed and (not (" << ConvolutionalLayer::getOutputFeatureCounterOutputSignalName() << "_Equals_Max and " << ConvolutionalLayer::getInputFeatureCounterOutputSignalName() << "_Equals_Max)));" << endl;
    }

    void ConvolutionalLayer::createDMAControllerSerial(Target* target)
    {
        unsigned int weightsPerAccess = this->myArguments->getCoreSize() * this->myArguments->getCoreSize();
        // before: counter width
        WeightFetcher* fetch = new WeightFetcher(target,weightsPerAccess,this->myArguments->getWeightWordSize(),this->myArguments->getStartAddress(),this->myArguments->getInputDepth()-1,this->myArguments->getNumberOfOutputFeatures()-1);
        // add component declaration
        addSubComponent(fetch);
        // connection with global DMA controller
        inPortMap(fetch,WeightFetcher::getDataInPortName(),Layer::getWeightInputName());
        inPortMap(fetch,WeightFetcher::getValidInPortName(),Layer::getWeightsValidInputName());
        inPortMap(fetch,WeightFetcher::getLastInPortName(),Layer::getLastWeightsInputName());
        outPortMap(fetch,WeightFetcher::getNewReadAccessPortName(),Layer::getNewReadAccessOutputName(),false);
        outPortMap(fetch,WeightFetcher::getNumberOfBytesPortName(),Layer::getNumberOfBytesOutputName(),false);
        outPortMap(fetch,WeightFetcher::getStartAddressPortName(),Layer::getNextStartAddressOutputName(),false);
        // connection with rest of the layer
        if(this->myArguments->getInputDepth()>1)
        {
            inPortMap(fetch,WeightFetcher::getInnerCounterPortName(),ConvolutionalLayer::getInputFeatureCounterOutputSignalName());
        }
        if(this->myArguments->getNumberOfOutputFeatures()>1)
        {
            inPortMap(fetch,WeightFetcher::getOuterCounterPortName(),ConvolutionalLayer::getOutputFeatureCounterOutputSignalName());
        }
        inPortMap(fetch,WeightFetcher::getStartNewReadPortName(),"Read_New_Weights");
        inPortMap(fetch,WeightFetcher::getWeightsArrivedAtConvCorePortName(),"Conv_Core_Weight_Enable");
        outPortMap(fetch,WeightFetcher::getNextWeightsValidPortName(),"New_Weights_Ready",false);
        for(unsigned int i=0; i<weightsPerAccess; i++)
        {
            outPortMap(fetch,fetch->getNextWeightsPortName(i),ConvolutionalLayer::getConvCoreWeightSignalName(i),false);
        }

        // create instance in vhdl code
        this->vhdl << instance(fetch,"WeightFetcher");
    }

    void ConvolutionalLayer::createCountersSerial(Target* target)
    {
        if(this->myArguments->getInputDepth() > 1 || this->myArguments->getNumberOfOutputFeatures() > 1)
        {
            // enable for modulo counters on rising edge of finished_delayed-signal
            Register* finishedReg = new Register(target,1,1);
            addSubComponent(finishedReg);
            inPortMap(finishedReg,"X","finished_delayed");
            outPortMap(finishedReg,"R","finished_delayed_reg",true);
            this->vhdl << instance(finishedReg,"Finished_Delayed_Reg_instance");
            this->vhdl << "finished_delayed_fall <= (not finished_delayed) and finished_delayed_reg;" << endl;
            this->vhdl << "counter_reset <= \"0\";" << endl;
        }

        if(this->myArguments->getInputDepth() > 1)
        {
            // input feature counter
            ModuloCounter* inC = new ModuloCounter(target,this->myArguments->getInputDepth());
            addSubComponent(inC);
            inPortMap(inC,"enable","finished_delayed_fall");
            inPortMap(inC,"manualReset","counter_reset");
            outPortMap(inC,"counter",ConvolutionalLayer::getInputFeatureCounterOutputSignalName(),false);
            this->vhdl << instance(inC,"Input_Feature_Counter");


            this->vhdl << "process(" << ConvolutionalLayer::getInputFeatureCounterOutputSignalName() << ")" << endl;
            this->vhdl << "begin" << endl;
            this->vhdl << tab << "if(unsigned(" << ConvolutionalLayer::getInputFeatureCounterOutputSignalName() << ") = " << this->myArguments->getInputDepth()-1 << ") then" << endl;
            this->vhdl << tab << tab << ConvolutionalLayer::getInputFeatureCounterOutputSignalName() << "_Equals_Max <= \"1\";" << endl;
            this->vhdl << tab << "else" << endl;
            this->vhdl << tab << tab << ConvolutionalLayer::getInputFeatureCounterOutputSignalName() << "_Equals_Max <= \"0\";" << endl;
            this->vhdl << tab << "end if;" << endl;
            this->vhdl << "end process;" << endl;

            if(this->inputMemoryParallelAccess==false)
            {
                this->vhdl << Layer::getGetNewDataName(0) << "_reset <= finished_delayed_fall;" << endl;
                // create lut that selects the correct input address based on the input feature counter
                this->createGetNewDataAddressLut(target);
            }
        }
        else
        {
            this->vhdl << ConvolutionalLayer::getInputFeatureCounterOutputSignalName() << " <= \"0\";" << endl;
            this->vhdl << ConvolutionalLayer::getInputFeatureCounterOutputSignalName() << "_Equals_Max" << " <= \"1\";" << endl;
        }

        if(this->myArguments->getNumberOfOutputFeatures() > 1)
        {
            // create enable-signal for output feature counter
            // ConvolutionalLayer::getOutputFeatureCounterOutputSignalName()+"_Equals_Max"
            this->vhdl << "process(" << ConvolutionalLayer::getOutputFeatureCounterOutputSignalName() << ")" << endl;
            this->vhdl << "begin" << endl;
            this->vhdl << tab << "if(unsigned(" << ConvolutionalLayer::getOutputFeatureCounterOutputSignalName() << ") = " << this->myArguments->getNumberOfOutputFeatures()-1 << ") then" << endl;
            this->vhdl << tab << tab << ConvolutionalLayer::getOutputFeatureCounterOutputSignalName() << "_Equals_Max <= \"1\";" << endl;
            this->vhdl << tab << "else" << endl;
            this->vhdl << tab << tab << ConvolutionalLayer::getOutputFeatureCounterOutputSignalName() << "_Equals_Max <= \"0\";" << endl;
            this->vhdl << tab << "end if;" << endl;
            this->vhdl << "end process;" << endl;

            this->vhdl << "Enable_Output_Feature_Counter <= finished_delayed_fall and " << ConvolutionalLayer::getInputFeatureCounterOutputSignalName() << "_Equals_Max;" << endl;
            // also set getNewData_reset
            if(this->inputMemoryParallelAccess==true)
            {
                for(int i=0; i<this->myArguments->getInputDepth(); i++)
                {
                    this->vhdl << Layer::getGetNewDataName(i) << "_reset <= Enable_Output_Feature_Counter;" << endl;
                }
            }

            // output feature counter
            ModuloCounter* outC = new ModuloCounter(target,this->myArguments->getNumberOfOutputFeatures());
            addSubComponent(outC);
            inPortMap(outC,"enable","Enable_Output_Feature_Counter");
            inPortMap(outC,"manualReset","counter_reset");
            outPortMap(outC,"counter",ConvolutionalLayer::getOutputFeatureCounterOutputSignalName(),false);
            this->vhdl << instance(outC,"Output_Feature_Counter");
        }
        else
        {
            this->vhdl << ConvolutionalLayer::getOutputFeatureCounterOutputSignalName() << " <= \"0\";" << endl;
            this->vhdl << ConvolutionalLayer::getOutputFeatureCounterOutputSignalName()+"_Equals_Max" << " <= \"1\";" << endl;
            // also set getNewData_reset
            if(this->inputMemoryParallelAccess==true)
            {
                for(int i=0; i<this->myArguments->getInputDepth(); i++)
                {
                    this->vhdl << Layer::getGetNewDataName(i) << "_reset <= \"0\";" << endl;
                }
            }
            else
            {
                this->vhdl << Layer::getGetNewDataName(0) << "_reset <= \"0\";" << endl;
            }

        }
    }

    void ConvolutionalLayer::createActivationFunctionSerial(Target* target)
    {
        // round
        string signalNameFrom = "Bitheap_out";
        string signalNameTo = "Bitheap_out_rounded";
        int wordSizeFrom = getSignalByName("Bitheap_out")->width();
        int fractionFrom = this->myArguments->getFraction() + this->myArguments->getWeightFraction();
        this->roundOutput(target,wordSizeFrom,fractionFrom,signalNameFrom,signalNameTo,"Rounding_instance",this->roundingType);

        // create relu and mux if needed
        if(this->myArguments->getActivationFunction()=="relu")
        {
            ReLU* relu = new ReLU(target,this->myArguments->getWordSize());
            addSubComponent(relu);
            inPortMap(relu,"X",signalNameTo);
            outPortMap(relu,"R",declare("Relu_out",this->myArguments->getWordSize()),false);
            this->vhdl << instance(relu,"ReLU");

            GenericMux* mux = new GenericMux(target,this->myArguments->getWordSize(),2);
            addSubComponent(mux);
            inPortMap(mux,mux->getSelectName(),ConvolutionalLayer::getInputFeatureCounterOutputSignalName()+"_Equals_Max");
            inPortMap(mux,mux->getInputName(0),"Bitheap_out_rounded");
            inPortMap(mux,mux->getInputName(1),"Relu_out");
            outPortMap(mux,mux->getOutputName(),"R_temp",false);
            this->vhdl << instance(mux,"ActivationFunctionMux");
        }
        else
        {
            this->vhdl << "R_temp <= Bitheap_out_rounded;" << endl;
        }

        // connect {R_0, ... , R_N} with R_temp
        if(this->outputMemoryParallelAccess==true)
        {
            for(int i=0; i<this->myArguments->getNumberOfOutputFeatures(); i++)
            {
                this->vhdl << this->getOutputSignalName(i) << " <= R_temp;" << endl;
            }
        }
        else
        {
            this->vhdl << this->getOutputSignalName(0) << " <= R_temp;" << endl;
        }

    }

    void ConvolutionalLayer::createLutForIntermediateReset(Target* target)
    {
        if(this->addressWidth!=0) // only build LUT if it's needed
        {
            map<unsigned int, unsigned int> lutData;
            unsigned int address = 0;
            for (int i = 0; i < this->myArguments->getNumberOfOutputFeatures(); i++) {
                lutData[i] = address;
                address += (this->outputWidth * this->outputHeight);
            }
            GenericLut *lut = new GenericLut(target, "Reset_Lut" + this->myArguments->getId(), lutData,
                                             this->getOutputCounterWidth(), this->addressWidth);
            addSubComponent(lut);
            for (unsigned int i = 0; i < this->getOutputCounterWidth(); i++) {
                this->vhdl << this->declare("Reset_Lut_i" + to_string(i), 1, false) << " <= "
                           << ConvolutionalLayer::getOutputFeatureCounterOutputSignalName() << "(" << i << ");" << endl;
                inPortMap(lut, "i" + to_string(i), "Reset_Lut_i" + to_string(i));
            }
            for (unsigned int i = 0; i < this->addressWidth; i++) {
                outPortMap(lut, "o" + to_string(i), "Reset_Lut_o" + to_string(i), true);
                this->vhdl << "reset_address(" << i << ") <= Reset_Lut_o" << i << ";" << endl;
            }
            this->vhdl << instance(lut, "Reset_Lut_instance");
        }
    }

    int ConvolutionalLayer::getInputCounterWidth() const
    {
        if(this->myArguments->getInputDepth()==1) return 1;
        return ceil(log2(this->myArguments->getInputDepth()));
    }

    int ConvolutionalLayer::getOutputCounterWidth() const
    {
        if(this->myArguments->getNumberOfOutputFeatures()==1) return 1;
        return ceil(log2(this->myArguments->getNumberOfOutputFeatures()));
    }


    void ConvolutionalLayer::generateVHDLCode(Target* target)
    {
        if(myArguments->getPaddingLeft()<0)
        {
            stringstream e;
            e << "Padding < 0 is not supported!";
            THROWERROR(e.str());
        }
        if(myArguments->getPaddingRight()<0)
        {
            this->myArguments->setPaddingRight(myArguments->getPaddingLeft());
        }
        if(myArguments->getPaddingTop()<0)
        {
            this->myArguments->setPaddingTop(myArguments->getPaddingLeft());
        }
        if(myArguments->getPaddingBot()<0)
        {
            this->myArguments->setPaddingBot(myArguments->getPaddingLeft());
        }


        this->useNumericStd();

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="ConvolutionalLayer";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // definition of the name of the operator
        ostringstream name;
        name << "ConvolutionalLayer_" << myArguments->getId();
        setName(name.str());

        // add a signal to start a new step
        addInput("newStep",1);
        addOutput("finished",1);

        if(this->myArguments->getCalcAllParallel()==true)
        {
            this->buildParallelCalculation(target);
        }
        else
        {
            this->buildSerialCalculation(target);
        }

    }

    void ConvolutionalLayer::createGetNewDataAddressLut(Target *target)
    {
        // calculate reset addresses
        map<unsigned int, unsigned int> addressMap;
        int inputSize = this->myArguments->getInputHeight() * this->myArguments->getInputWidth();
        for(unsigned int i=0; i<this->myArguments->getInputDepth(); i++)
        {
            if(i<(this->myArguments->getInputDepth()-1))
            {
                // set the address on the first address of the next feature
                addressMap[i] = (i+1) * inputSize;
            }
            else
            {
                // the next calculation will be the first input feature of the next output feature => set the address back to zero
                addressMap[i] = 0;
            }
        }

        // create and connect lut
        int wOut_Lut = (ceil(log2(this->myArguments->getInputWidth()*this->myArguments->getInputHeight()*this->myArguments->getInputDepth())));
        GenericLut* lut = new GenericLut(target,"Get_New_Data_Reset_Address_Lut_"+this->myArguments->getId(),addressMap,this->getInputCounterWidth(),wOut_Lut);
        addSubComponent(lut);
        for(int i=0; i<this->getInputCounterWidth(); i++)
        {
            this->vhdl << declare(ConvolutionalLayer::getInputFeatureCounterOutputSignalName()+to_string(i)+"_std",1,false)
                       << " <= " << ConvolutionalLayer::getInputFeatureCounterOutputSignalName() << "(" << i << ");" << endl;
            inPortMap(lut,"i"+to_string(i),ConvolutionalLayer::getInputFeatureCounterOutputSignalName()+to_string(i)+"_std");
        }
        for(int i=0; i<wOut_Lut; i++)
        {
            this->vhdl << Layer::getGetNewDataName(0) << "_reset_address(" << i << ") <= " << declare("Get_New_Data_Reset_Address_Lut_out_"+to_string(i),1,false) << ";" << endl;
            outPortMap(lut,"o"+to_string(i),"Get_New_Data_Reset_Address_Lut_out_"+to_string(i),false);
        }
        this->vhdl << instance(lut,"Get_New_Data_Reset_Address_Lut_instance");
    }

}//namespace flopoco
