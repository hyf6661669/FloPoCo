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
#include "WindowShiftRegister.hpp"
#include "PaddingGenerator.hpp"
#include "ReLU.hpp"
#include "IntAddSubCmp/IntAdderTree.hpp"
#include "Register.hpp"


using namespace std;
namespace flopoco {




    ConvolutionalLayer::ConvolutionalLayer(Target* target, unsigned int wordSize_, unsigned int fraction_, unsigned int horizontalSize_, unsigned int verticalSize_, unsigned int windowSize_, unsigned int numberOfFeatures_, unsigned int numberOfInputFeatures_, vector <vector <vector <double> > > weights_, unsigned int weightWordSize_, unsigned int weightFraction_, int paddingLeft_, string paddingType_, unsigned int stride_, bool useAdderTree_, bool roundAfterConvCore_, string roundingType_, string activationFunction_, int paddingRight_, int paddingTop_, int paddingBot_, string id_) :
        Operator(target), wordSize(wordSize_), fraction(fraction_), horizontalSize(horizontalSize_), verticalSize(verticalSize_), windowSize(windowSize_), numberOfFeatures(numberOfFeatures_), numberOfInputFeatures(numberOfInputFeatures_), weights(weights_), weightWordSize(weightWordSize_), weightFraction(weightFraction_), paddingLeft(paddingLeft_), paddingRight(paddingRight_), paddingTop(paddingTop_), paddingBot(paddingBot_), paddingType(paddingType_), stride(stride_), useAdderTree(useAdderTree_), roundAfterConvCore(roundAfterConvCore_), roundingType(roundingType_), activationFunction(activationFunction_), id(id_) {

        if(paddingLeft<0)
        {
            stringstream e;
            e << "Padding < 0 is not supported!";
            THROWERROR(e);
        }
        if(paddingRight<0)
        {
            paddingRight = paddingLeft;
        }
        if(paddingTop<0)
        {
            paddingTop = paddingLeft;
        }
        if(paddingBot<0)
        {
            paddingBot = paddingLeft;
        }

        this->useNumericStd();

		// definition of the source file name, used for info and error reporting using REPORT 
        srcFileName="ConvolutionalLayer";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

		// definition of the name of the operator
		ostringstream name;
        name << "ConvolutionalLayer_" << id;
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

        for(unsigned int featureCounter=0; featureCounter<numberOfFeatures; featureCounter++)
        {
            vhdl << "validData_o_" << featureCounter << " <= validData_temp;" << endl;
        }

        // activation function (atm only ReLU is supported)
        cout << "###   building activation function \n";
        this->buildActivationFunction(target);

    }

    void ConvolutionalLayer::buildShiftAndPad(Target* target)
    {

        cout << "###   ConvLayer.buildShiftAndPad \n";
        for(unsigned int featureCounter=0; featureCounter<numberOfInputFeatures; featureCounter++)
        {
            cout << "###       featureCounter=" << featureCounter << "\n";

            cout << "###       declaring inputs \n";
            // declare inputs
            addInput("X"+to_string(featureCounter), wordSize);
            addInput("validData_i_"+to_string(featureCounter),1);
            addOutput("getNewData"+to_string(featureCounter),1);

            // shift register
            cout << "###       building shift register \n";
            WindowShiftRegister* shiftReg = new WindowShiftRegister(target,wordSize,windowSize,horizontalSize);
            addSubComponent(shiftReg);
            inPortMap(shiftReg,"X","X"+to_string(featureCounter));
            inPortMap(shiftReg,"enable","validData_i_"+to_string(featureCounter));
            inPortMap(shiftReg,"newStep","newStep");
            for(unsigned int portC=0;portC<windowSize*windowSize;portC++)
            {
                outPortMap(shiftReg,"R"+to_string(portC),"PaddingGenInput_inputFeature_"+to_string(featureCounter)+"_number_"+to_string(portC),true);
            }
            vhdl << instance(shiftReg,"WindowShiftRegister"+to_string(featureCounter)+"_instance") << endl;

            // padding generator

            cout << "###       building padding generator \n";
            PaddingGenerator* padGen = new PaddingGenerator(target,wordSize,windowSize,horizontalSize,verticalSize,paddingTop,stride,paddingType,paddingBot,paddingLeft,paddingRight,(featureCounter==0?true:false));
            addSubComponent(padGen);
            for(unsigned int portC=0; portC<windowSize*windowSize; portC++)
            {
                inPortMap(padGen,"X"+to_string(portC),"PaddingGenInput_inputFeature_"+to_string(featureCounter)+"_number_"+to_string(portC));
                outPortMap(padGen,"R"+to_string(portC),"PaddingGenOutput_inputFeature_"+to_string(featureCounter)+"_number_"+to_string(portC),true);
            }
            inPortMap(padGen,"validData_i","validData_i_"+to_string(featureCounter));
            inPortMap(padGen,"newStep","newStep");
            outPortMap(padGen,"getNewData","getNewData"+to_string(featureCounter),false);
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
        unsigned int convCoreIdCounter=0;
        for(unsigned int featureCounter=0; featureCounter<numberOfInputFeatures; featureCounter++)
        {
            cout << "###       featureCounter=" << featureCounter << "\n";
            // one ConvCore for each OUTPUT feature
            // (get pipelinestages from all ConvCores for each OUTPUT-feature and insert additional Delays if needed)
            vector <ConvolutionalCoreSimple*> tempConvCores;
            this->ConvCores.push_back(tempConvCores);
            unsigned int tempOutWordSize;
            unsigned int tempOutFraction;
            for(unsigned int outFeatureCounter=0; outFeatureCounter<numberOfFeatures; outFeatureCounter++)
            {
                cout << "###       outFeatureCounter=" << outFeatureCounter << "\n";
                cout << "###       ConvCore \n";
                cout << "###           wordSize=" << wordSize << endl;
                cout << "###           fraction=" << fraction << endl;
                cout << "###           weightWordSize=" << weightWordSize << endl;
                cout << "###           weightFraction=" << weightFraction << endl;
                cout << "###           windowSize=" << windowSize << endl;
                cout << "###           useAdderTree=" << useAdderTree << endl;
                cout << "###           roundAfterConvCore=" << roundAfterConvCore << endl;
                cout << "###           roundingType=" << roundingType << endl;
                cout << "###           convCoreIdCounter=" << convCoreIdCounter << endl;
                cout << "###           weights:" << endl;
                for(auto it : weights[featureCounter][outFeatureCounter])
                {
                    cout << "###               " << it << endl;
                }
                ConvolutionalCoreSimple* convCore = new ConvolutionalCoreSimple(target,wordSize,fraction,weightWordSize,weightFraction,windowSize,weights[featureCounter][outFeatureCounter],useAdderTree,(roundAfterConvCore==true?roundingType:""),to_string(convCoreIdCounter));
                cout << "###       after ConvCore \n";
                convCoreIdCounter++;
                this->addSubComponent(convCore);
                for(unsigned int portC=0; portC<windowSize*windowSize; portC++)
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

    void ConvolutionalLayer::buildAdderTrees(Target* target)
    {
        //gen adder trees
        for(unsigned int outFeatureCounter=0; outFeatureCounter<this->numberOfFeatures; outFeatureCounter++)
        {
            cout << "###       outFeatureCounter=" << outFeatureCounter << "\n";
            //validData_o
            addOutput("validData_o_"+to_string(outFeatureCounter),1);

            //one adder tree for each output feature
            cout << "###       IntAdderTree \n";
            IntAdderTree* intAdd = new IntAdderTree(target,this->ConvCores[0][outFeatureCounter]->getOutputWordSize(),this->numberOfInputFeatures, "bitheap",false);
            this->AdderTrees.push_back(intAdd);
            addSubComponent(intAdd);
            for(unsigned int portC=0; portC<numberOfInputFeatures; portC++)
            {
                inPortMap(intAdd,"X"+to_string(portC+1),"ConvCoreOutput_inputFeature_"+to_string(portC)+"_outputFeature_"+to_string(outFeatureCounter));
            }
            outPortMap(intAdd,"Y","AdderTreeOutputFeature_"+to_string(outFeatureCounter)+"_out",true);
            vhdl << instance(intAdd,"AdderTreeOutputFeature_"+to_string(outFeatureCounter)+"_instance");

            // round back to original word size
        }

        //sync
        unsigned int outFCCounterAdd = 0;
        for(unsigned int i=0; i<this->AdderTrees.size(); i++)
        {
            if(this->AdderTrees[i]->getPipelineDepth()>=AdderTrees[outFCCounterAdd]->getPipelineDepth())
            {
                outFCCounterAdd = i;
            }
        }
        string tempSignalName="AdderTreeOutputFeature_"+to_string(outFCCounterAdd)+"_out"; // this is the output-signal from the AdderTree with the most internal pipelinestages
        this->syncCycleFromSignal(tempSignalName);

        //round
        for(unsigned int outFeatureCounter=0; outFeatureCounter<this->numberOfFeatures; outFeatureCounter++)
        {
            cout << "###       Rounding back to original word size \n";
            roundOutput(this->getSignalByName("AdderTreeOutputFeature_"+to_string(outFeatureCounter)+"_out")->width(),this->ConvCores[0][outFeatureCounter]->getOutputFraction(),"AdderTreeOutputFeature_"+to_string(outFeatureCounter)+"_out","feature"+to_string(outFeatureCounter)+"_rounded",this->roundingType);
        }
    }

    void ConvolutionalLayer::roundOutput(unsigned int wordSizeFrom, unsigned int fractionFrom, string signalNameFrom, string signalNameTo, string round)
    {
        cout << "###   ConvLayer.roundOutput \n";
        if(round=="Truncation")
        {
            int MSBFrom = wordSizeFrom - fractionFrom;
            int MSBTo = this->wordSize-this->fraction;
            int from = wordSizeFrom-(MSBFrom-MSBTo)-1;
            int downto = fractionFrom-this->fraction;
            declare(signalNameTo,wordSize);

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
        if(this->activationFunction=="ReLU")
        {
            for(unsigned int outFeatureCounter=0; outFeatureCounter<this->numberOfFeatures; outFeatureCounter++)
            {
                ReLU* relu = new ReLU(target, this->wordSize);
                addSubComponent(relu);
                inPortMap(relu,"X","feature"+to_string(outFeatureCounter)+"_rounded");
                outPortMap(relu,"R","ReLU"+to_string(outFeatureCounter)+"_out",true);
                vhdl << instance(relu,"ReLU"+to_string(outFeatureCounter)+"_instance");

                addOutput("R"+to_string(outFeatureCounter),this->wordSize);
                vhdl << "R" << outFeatureCounter << " <= ReLU" << outFeatureCounter << "_out;" << endl;
            }
        }
        else
        {
            for(unsigned int outFeatureCounter=0; outFeatureCounter<this->numberOfFeatures; outFeatureCounter++)
            {
                addOutput("R"+to_string(outFeatureCounter),this->wordSize);
                vhdl << "R" << outFeatureCounter << " <= feature" << outFeatureCounter << "_rounded;" << endl;
            }
        }
    }

}//namespace flopoco
