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

// include the header of this Operator and the header of all nother ecessary Operators
#include "ConvolutionalLayer.hpp"
#include "ConvolutionalCoreSimple.hpp"
#include "NeuralNetworks/Utility/WindowShiftRegister.hpp"
#include "NeuralNetworks/Utility/PaddingGenerator.hpp"
#include "NeuralNetworks/ActivationFunctions/ReLU.hpp"
//#include "IntAddSubCmp/IntAdderTree.hpp"
#include "NeuralNetworks/Utility/Register.hpp"


using namespace std;
namespace flopoco {




    ConvolutionalLayer::ConvolutionalLayer(Target* target, int wordSize_, int fraction_, int horizontalSize_, int verticalSize_, int windowSize_, int numberOfOutputFeatures_, int numberOfInputFeatures_, vector <vector <vector <double> > > weights_, int weightWordSize_, int weightFraction_, int paddingLeft_, string paddingType_, bool inputFeaturesParallel_, bool outputFeaturesParallel_, string id_, int stride_, bool useAdderTree_, bool roundAfterConvCore_, string roundingType_, string activationFunction_, int paddingRight_, int paddingTop_, int paddingBot_) :
        Layer(target), useAdderTree(useAdderTree_), roundAfterConvCore(roundAfterConvCore_), roundingType(roundingType_) {

//        this->wordSize=wordSize_;
//        this->fraction=fraction_;
//        this->horizontalSize=horizontalSize_;
//        this->verticalSize=verticalSize_;
//        this->numberOfOutputFeatures=numberOfOutputFeatures_;
//        this->numberOfInputFeatures=numberOfInputFeatures_;
//        this->weightWordSize=weightWordSize_;
//        this->weightFraction=weightFraction_;
//        this->inputFeaturesParallel=inputFeaturesParallel_;
//        this->outputFeaturesParallel=outputFeaturesParallel_;
//        this->activationFunction=activationFunction_;
//        this->layerType="Convolutional";
//        this->myArguments = new LayerArguments();
        this->myArguments->setWordSize(wordSize_); //
        this->myArguments->setFraction(fraction_); //
        this->myArguments->setInputWidth(horizontalSize_); //
        this->myArguments->setInputHeight(verticalSize_); //
        this->myArguments->setNumberOfOutputFeatures(numberOfOutputFeatures_); //
        this->myArguments->setInputDepth(numberOfInputFeatures_); //
        this->myArguments->setWeightWordSize(weightWordSize_); //
        this->myArguments->setWeightFraction(weightFraction_); //
        this->myArguments->setInputFeaturesParallel(inputFeaturesParallel_); //
        this->myArguments->setOutputFeaturesParallel(outputFeaturesParallel_); //
        this->myArguments->setActivationFunction(activationFunction_); //
        this->myArguments->setLayerType("Convolutional"); //
        this->myArguments->setCoreSize(windowSize_); //
        this->myArguments->setConvWeights(weights_); //
        this->myArguments->setPaddingLeft(paddingLeft_); //
        this->myArguments->setPaddingRight(paddingRight_); //
        this->myArguments->setPaddingTop(paddingTop_); //
        this->myArguments->setPaddingBot(paddingBot_); //
        this->myArguments->setPaddingType(paddingType_); //
        this->myArguments->setStride(stride_); //
        this->myArguments->setId(id_); //

        generateVHDLCode(target);
    }

    ConvolutionalLayer::ConvolutionalLayer(Target *target, LayerArguments *la, bool useAdderTree_, bool roundAfterConvCore_, string roundingType_) :
        Layer(target, la), useAdderTree(useAdderTree_), roundAfterConvCore(roundAfterConvCore_), roundingType(roundingType_)
    {
        generateVHDLCode(target);
    }

    string ConvolutionalLayer::getOutputSignalName(int feature)
    {
        if(this->myArguments->getNumberOfOutputFeatures()<=feature)
        {
            cout << "ConvolutionalLayer.getOutputSignalName: requested output signal doesn't exist" << endl;
        }
        return "R"+to_string(feature);
    }

    string ConvolutionalLayer::getInputSignalName(int feature)
    {
        if(this->myArguments->getInputDepth()<=feature)
        {
            cout << "ConvolutionalLayer.getInputSignalName: requested input signal doesn't exist" << endl;
        }
        return "X"+to_string(feature);
    }

    void ConvolutionalLayer::buildShiftAndPad(Target* target)
    {

        cout << "###   ConvLayer.buildShiftAndPad \n";
        for(int featureCounter=0; featureCounter<myArguments->getInputDepth(); featureCounter++)
        {
            cout << "###       featureCounter=" << featureCounter << "\n";

            cout << "###       declaring inputs \n";
            // declare inputs
            addInput("X"+to_string(featureCounter), myArguments->getWordSize());
            addInput("validData_i_"+to_string(featureCounter),1);
            addOutput("getNewData_"+to_string(featureCounter),1);

            // shift register
            cout << "###       building shift register \n";
            WindowShiftRegister* shiftReg = new WindowShiftRegister(target,myArguments->getWordSize(),myArguments->getCoreSize(),myArguments->getInputWidth());
            addSubComponent(shiftReg);
            inPortMap(shiftReg,"X","X"+to_string(featureCounter));
            inPortMap(shiftReg,"enable","validData_i_"+to_string(featureCounter));
            inPortMap(shiftReg,"newStep","newStep");
            for(int portC=0;portC<myArguments->getCoreSize()*myArguments->getCoreSize();portC++)
            {
                outPortMap(shiftReg,"R"+to_string(portC),"PaddingGenInput_inputFeature_"+to_string(featureCounter)+"_number_"+to_string(portC),true);
            }
            vhdl << instance(shiftReg,"WindowShiftRegister"+to_string(featureCounter)+"_instance") << endl;

            // padding generator

            cout << "###       building padding generator \n";
            PaddingGenerator* padGen = new PaddingGenerator(target,myArguments->getWordSize(),myArguments->getCoreSize(),myArguments->getInputWidth(),myArguments->getInputHeight(),myArguments->getPaddingTop(),myArguments->getStride(),myArguments->getPaddingType(),myArguments->getPaddingBot(),myArguments->getPaddingLeft(),myArguments->getPaddingRight(),(featureCounter==0?true:false));
            addSubComponent(padGen);
            for(int portC=0; portC<myArguments->getCoreSize()*myArguments->getCoreSize(); portC++)
            {
                inPortMap(padGen,"X"+to_string(portC),"PaddingGenInput_inputFeature_"+to_string(featureCounter)+"_number_"+to_string(portC));
                outPortMap(padGen,"R"+to_string(portC),"PaddingGenOutput_inputFeature_"+to_string(featureCounter)+"_number_"+to_string(portC),true);
            }
            inPortMap(padGen,"validData_i","validData_i_"+to_string(featureCounter));
            inPortMap(padGen,"newStep","newStep");
            outPortMap(padGen,"getNewData","getNewData_"+to_string(featureCounter),false);
            if(featureCounter==0)
            {
                outPortMap(padGen,"validData_o","validData_temp",true);
                outPortMap(padGen,"finished","finished_tmp",true);
            }
            vhdl << instance(padGen,"PaddingGen_"+to_string(featureCounter)+"_instance");

        }
    }

    void ConvolutionalLayer::buildConvCores(Target* target)
    {
        cout << "###   ConvolutionalLayer.buildConvCores \n";
        int convCoreIdCounter=0;
        for(int featureCounter=0; featureCounter<myArguments->getInputDepth(); featureCounter++)
        {
            cout << "###       featureCounter=" << featureCounter << "\n";
            // one ConvCore for each OUTPUT feature
            // (get pipelinestages from all ConvCores for each OUTPUT-feature and insert additional Delays if needed)
            vector <ConvolutionalCoreSimple*> tempConvCores;
            this->ConvCores.push_back(tempConvCores);
            int tempOutWordSize;
            int tempOutFraction;
            for(int outFeatureCounter=0; outFeatureCounter<myArguments->getNumberOfOutputFeatures(); outFeatureCounter++)
            {
                cout << "###       outFeatureCounter=" << outFeatureCounter << "\n";
                cout << "###       ConvCore \n";
                cout << "###           wordSize=" << myArguments->getWordSize() << endl;
                cout << "###           fraction=" << myArguments->getFraction() << endl;
                cout << "###           weightWordSize=" << myArguments->getWeightWordSize() << endl;
                cout << "###           weightFraction=" << myArguments->getWeightFraction() << endl;
                cout << "###           windowSize=" << myArguments->getCoreSize() << endl;
                cout << "###           useAdderTree=" << useAdderTree << endl;
                cout << "###           roundAfterConvCore=" << roundAfterConvCore << endl;
                cout << "###           roundingType=" << roundingType << endl;
                cout << "###           convCoreIdCounter=" << convCoreIdCounter << endl;
                cout << "###           weights:" << endl;
                for(auto it : myArguments->getConvWeights()[featureCounter][outFeatureCounter])
                {
                    cout << "###               " << it << endl;
                }
                ConvolutionalCoreSimple* convCore = new ConvolutionalCoreSimple(target,myArguments->getWordSize(),myArguments->getFraction(),myArguments->getWeightWordSize(),myArguments->getWeightFraction(),myArguments->getCoreSize(),myArguments->getConvWeights()[featureCounter][outFeatureCounter],useAdderTree,(roundAfterConvCore==true?roundingType:""),to_string(convCoreIdCounter));
                cout << "###       after ConvCore \n";
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
                    cout << "###           tempOutWordSize=" << tempOutWordSize << endl;
                    cout << "###           convCore->getOutputWordSize()=" << convCore->getOutputWordSize() << endl;
                    cout << "###           tempOutFraction=" << tempOutFraction << endl;
                    cout << "###           convCore->getOutputFraction()=" << convCore->getOutputFraction() << endl;
                    if(tempOutWordSize!=convCore->getOutputWordSize() || tempOutFraction!=convCore->getOutputFraction())
                    {
                        cout << "###       FEHLERMELDUNG!!!!! \n";
                        stringstream e;
                        e << "The Convolutional Cores have a different outputWordSize!";
                        THROWERROR(e);
                    }
                }
                tempOutWordSize=convCore->getOutputWordSize();
                tempOutFraction=convCore->getOutputFraction();
            }
        }

        // delay data and flag-signals for all convCores with pipeline stages < maxPipeline stages
        int inputFCCounter=0;
        int outputFCCounter=0;
        for(int inputFC=0; inputFC<this->ConvCores.size(); inputFC++)
        {
            for(int outputFC=0; outputFC<this->ConvCores[inputFC].size(); outputFC++)
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

    void ConvolutionalLayer::buildAdderTrees(Target* target)
    {
        //gen adder trees
        for(int outFeatureCounter=0; outFeatureCounter<this->myArguments->getNumberOfOutputFeatures(); outFeatureCounter++)
        {
            cout << "###       outFeatureCounter=" << outFeatureCounter << "\n";
            //validData_o
            addOutput("validData_o_"+to_string(outFeatureCounter),1);

            //one adder tree for each output feature
            cout << "###       IntAdderTree \n";
            IntAdderTree* intAdd = new IntAdderTree(target,this->ConvCores[0][outFeatureCounter]->getOutputWordSize(),this->myArguments->getInputDepth(), "bitheap",false);
            this->AdderTrees.push_back(intAdd);
            addSubComponent(intAdd);
            for(int portC=0; portC<myArguments->getInputDepth(); portC++)
            {
                inPortMap(intAdd,"X"+to_string(portC+1),"ConvCoreOutput_inputFeature_"+to_string(portC)+"_outputFeature_"+to_string(outFeatureCounter));
            }
            outPortMap(intAdd,"Y","AdderTreeOutputFeature_"+to_string(outFeatureCounter)+"_out",true);
            vhdl << instance(intAdd,"AdderTreeOutputFeature_"+to_string(outFeatureCounter)+"_instance");

            // round back to original word size
        }

        //sync
        int outFCCounterAdd = 0;
        for(int i=0; i<this->AdderTrees.size(); i++)
        {
            if(this->AdderTrees[i]->getPipelineDepth()>=AdderTrees[outFCCounterAdd]->getPipelineDepth())
            {
                outFCCounterAdd = i;
            }
        }
        string tempSignalName="AdderTreeOutputFeature_"+to_string(outFCCounterAdd)+"_out"; // this is the output-signal from the AdderTree with the most internal pipelinestages
        this->syncCycleFromSignal(tempSignalName);

        //round
        for(int outFeatureCounter=0; outFeatureCounter<this->myArguments->getNumberOfOutputFeatures(); outFeatureCounter++)
        {
            cout << "###       Rounding back to original word size \n";
            roundOutput(this->getSignalByName("AdderTreeOutputFeature_"+to_string(outFeatureCounter)+"_out")->width(),this->ConvCores[0][outFeatureCounter]->getOutputFraction(),"AdderTreeOutputFeature_"+to_string(outFeatureCounter)+"_out","feature"+to_string(outFeatureCounter)+"_rounded",this->roundingType);
        }
    }

    void ConvolutionalLayer::roundOutput(int wordSizeFrom, int fractionFrom, string signalNameFrom, string signalNameTo, string round)
    {
        cout << "###   ConvLayer.roundOutput \n";
        if(round=="Truncation")
        {
            int MSBFrom = wordSizeFrom - fractionFrom;
            int MSBTo = this->myArguments->getWordSize()-this->myArguments->getFraction();
            int from = wordSizeFrom-(MSBFrom-MSBTo)-1;
            int downto = fractionFrom-this->myArguments->getFraction();
            declare(signalNameTo,myArguments->getWordSize());

            if(downto>0) // correct rounding when cutting off LSBs
            {
                this->vhdl << tab << declare(signalNameTo+"_roundingNumber",downto) << " <= " << signalNameFrom << "(" << downto-1 << ")";
                if(downto>1)
                {
                    this->vhdl << " & (" << downto-2 << " downto 0 => '0')";
                }
                this->vhdl <<";" << endl;

                this->vhdl << tab << declare(signalNameTo+"_beforeRound",wordSizeFrom) << " <= std_logic_vector(signed(" << signalNameFrom
                                  << ")+" << "signed("<< signalNameTo+"_roundingNumber" << "));" << endl;

                this->vhdl << tab << signalNameTo << " <= "<< signalNameTo << "_beforeRound(" << from << " downto " << downto << ");" << endl;
            }
            else if(downto==0) // only cut off MSBs
            {
                this->vhdl << tab << signalNameTo << " <= " << signalNameFrom << "(" << from << " downto " << downto << ");" << endl;
            }
            else
            {
                stringstream e;
                e << "Error while rounding, requested fraction is bigger than actual fraction! This should normally not happen";
                THROWERROR(e);
            }



        }
        else
        {
            cout << "ConvolutionalLayer.roundOutput: atm only truncation is implemented. FloPoCo is rounding with this mode instead" << endl;
            this->roundOutput(wordSizeFrom, fractionFrom, signalNameFrom, signalNameTo, "Truncation");
        }
    }

    void ConvolutionalLayer::buildActivationFunction(Target* target)
    {
        if(this->myArguments->getActivationFunction()=="ReLU")
        {
            for(int outFeatureCounter=0; outFeatureCounter<this->myArguments->getNumberOfOutputFeatures(); outFeatureCounter++)
            {
                ReLU* relu = new ReLU(target, this->myArguments->getWordSize());
                addSubComponent(relu);
                inPortMap(relu,"X","feature"+to_string(outFeatureCounter)+"_rounded");
                outPortMap(relu,"R","ReLU"+to_string(outFeatureCounter)+"_out",true);
                vhdl << instance(relu,"ReLU"+to_string(outFeatureCounter)+"_instance");

                addOutput("R"+to_string(outFeatureCounter),this->myArguments->getWordSize());
                vhdl << "R" << outFeatureCounter << " <= ReLU" << outFeatureCounter << "_out;" << endl;
            }
        }
        else
        {
            for(int outFeatureCounter=0; outFeatureCounter<this->myArguments->getNumberOfOutputFeatures(); outFeatureCounter++)
            {
                addOutput("R"+to_string(outFeatureCounter),this->myArguments->getWordSize());
                vhdl << "R" << outFeatureCounter << " <= feature" << outFeatureCounter << "_rounded;" << endl;
            }
        }
    }

    void ConvolutionalLayer::generateVHDLCode(Target* target)
    {
        if(myArguments->getPaddingLeft()<0)
        {
            stringstream e;
            e << "Padding < 0 is not supported!";
            THROWERROR(e);
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
        cout << "###   Adding Inputs (newStep and finished) \n";
        addInput("newStep",1);
        addOutput("finished",1);

        // build one convolutional chain for each INPUT feature!
        cout << "###   Building shift and pad \n";
        buildShiftAndPad(target);

        cout << "###   Building ConvCores \n";
        buildConvCores(target);

        cout << "###   FINISHED building ConvCores! \n";

        // Add all ConvCore-outputs for each OUTPUT-feature up and round back to original word size (while doing that, also create validData_o-signals because we need them anyways)
        cout << "###   Building Adder Trees for each outFeature \n";
        this->buildAdderTrees(target);
        cout << "###   FINISHED building Adder Trees \n";

        // delay finished-signal and validData-signal for x cycles, with x=maxNumberOfPipelineStages(all AdderTrees) and assign outputs
        cout << "###   building flag signals \n";
        vhdl << "finished <= finished_tmp;" << endl;

        for(int featureCounter=0; featureCounter<myArguments->getNumberOfOutputFeatures(); featureCounter++)
        {
            vhdl << "validData_o_" << featureCounter << " <= validData_temp;" << endl;
        }

        // activation function (atm only ReLU is supported)
        cout << "###   building activation function \n";
        this->buildActivationFunction(target);
    }

}//namespace flopoco
