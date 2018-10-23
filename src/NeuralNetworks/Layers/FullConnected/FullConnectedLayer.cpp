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

// include the header of the Operator and a few others
#include "FullConnectedLayer.hpp"
#include "FullConnectedCore.hpp"
#include "FullConnectedLayerController.hpp"

#include "NeuralNetworks/MemoryManagement/DmaAccess/WeightFetcher.hpp"
#include "NeuralNetworks/Utility/Register.hpp"
#include "NeuralNetworks/Utility/Rounding.hpp"
#include "NeuralNetworks/Utility/BitheapWrapper.hpp"
#include "NeuralNetworks/ActivationFunctions/ReLU.hpp"

#include "PrimitiveComponents/GenericLut.hpp"
#include "PrimitiveComponents/GenericMux.hpp"

using namespace std;
namespace flopoco {




    FullConnectedLayer::FullConnectedLayer(Target* target, NeuralNetwork* parent, LayerArguments* lA, unsigned int weightsPerDMAAccess_, char roundingType_) :
            Layer(target,lA,parent), weightsPerDMAAccess(weightsPerDMAAccess_), roundingType(roundingType_), fullConnectedCoreAdderPipelineDepth(-1) {

		// definition of the source file name, used for info and error reporting using REPORT 
		srcFileName="FullConnectedLayer";

        // author
        setCopyrightString("Nicolai Fiege, 2018");

		// definition of the name of the operator
		ostringstream name;
        name << "FullConnectedLayer_" << this->myArguments->getId();
        setName(name.str());

        // set input and output memory access to 'serial'
        this->inputMemoryParallelAccess=false;
        this->outputMemoryParallelAccess=false;

        this->generateVHDLCode(target);
    }

    void FullConnectedLayer::generateVHDLCode(Target* target)
    {
        // In/Out
        this->buildInputsAndOutputs(target);
        // declare signals
        this->declareSignals(target);
        // Bias LUT
        this->buildBiasLut(target);
        // build FullConnectedCore
        this->buildFullConnectedCore(target);
        // Controller
        this->buildController(target);
        // WeightFetcher
        this->buildWeightFetcher(target);
        // Weight Registers
        this->buildWeightRegisters(target);
        // Round
        this->roundOutput(target);
        // Activation Function
        this->buildActivationFunction(target);
    }

    string FullConnectedLayer::getInputSignalName(int number)
    {
        if(number>0) THROWERROR("Requested Input (" << number << ") doesn't exist!");
        return "X_0";
    }

    string FullConnectedLayer::getOutputSignalName(int number)
    {
        if(number>0) THROWERROR("Requested Output (" << number << ") doesn't exist!");
        return "R_0";
    }

    void FullConnectedLayer::buildInputsAndOutputs(Target* target)
    {
        //////////////////////////////////////////
        // Communication with Global Controller //
        //////////////////////////////////////////
        addInput("newStep",1);
        addOutput("finished",1);

        /////////////////////
        // DATA INPUT SIDE //
        /////////////////////
        // X
        addInput(this->getInputSignalName(0),this->myArguments->getWordSize());
        // validData_i
        addInput(Layer::getValidDataName(0)+"_i",1,true);
        // getNewData
        addOutput(Layer::getGetNewDataName(0),1,true);
        // getNewData_reset
        addOutput(Layer::getGetNewDataName(0)+"_reset",1,true);

        //////////////////////
        // DATA OUTPUT SIDE //
        //////////////////////
        // R
        addOutput(this->getOutputSignalName(0),this->myArguments->getWordSize());
        // validData_o
        addOutput(Layer::getValidDataName(0)+"_o",1,true);

        ///////////////////////////////////////////////////////
        // connect Weight fetcher with global DMA controller //
        ///////////////////////////////////////////////////////
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

    void FullConnectedLayer::declareSignals(Target* target)
    {
        this->weightsPerNeuron = this->myArguments->getInputHeight() * this->myArguments->getInputWidth() * this->myArguments->getInputDepth();
        this->neuronCounterWidth = ceil(log2(this->myArguments->getNumberOfOutputFeatures()));
        if(this->neuronCounterWidth==0) this->neuronCounterWidth=1;
        this->weightCounterWidth = ceil(log2(this->weightsPerNeuron));
        if(this->weightCounterWidth==0) this->weightCounterWidth=1;
        this->intermediateResultWidth = this->myArguments->getWordSize() + this->myArguments->getWeightWordSize() + weightCounterWidth;
        declare("Bias",this->intermediateResultWidth);
        declare("Neuron_Counter",this->neuronCounterWidth);
        declare("Weight_Counter",this->weightCounterWidth);
        declare("Intermediate_Result",this->intermediateResultWidth);
        declare("Intermediate_Result_Reg_Input",this->intermediateResultWidth);
        declare("Full_Connected_Core_Second_Input",this->intermediateResultWidth);
        declare("Full_Connected_Core_Weight_Input",this->myArguments->getWeightWordSize());
        declare("Full_Connected_Core_Output",this->intermediateResultWidth);
        declare("Weights_Are_Valid",1);
        declare("finished_temp",1);
        declare("Intermediate_Result_Reg_enable",1);
        declare("Intermediate_Result_Reg_reset",1);
        declare("Intermediate_Result_Reg_enable_reg",1);
        declare("Weight_Register_Enable",1);
        declare("Start_Fetching_New_Weights",1);
        declare("R_rounded",this->myArguments->getWordSize());
    }

    void FullConnectedLayer::buildBiasLut(Target* target)
    {
        // copy of bias vector
        vector<double> biases = this->myArguments->getBiases();
        declare("Bias_temp",this->myArguments->getWeightWordSize());
        if(biases.size()==0) THROWERROR("Bias Vector can't be empty!");
        if(biases.size()==1)
        {
            double scaledBias = biases[0] * pow(2, this->myArguments->getWeightFraction());
            this->vhdl << "Bias_temp <= std_logic_vector(signed(" << scaledBias << "," << this->myArguments->getWeightWordSize() << "));" << endl;
        }
        else
        {
            map<unsigned int, unsigned int> biasMap;
            for (unsigned int i = 0; i < biases.size(); i++) {
                double scaledBias = biases[i] * pow(2, this->myArguments->getWeightFraction());
                if ((double) abs(scaledBias) >= pow(2, 32)) {
                    THROWERROR("Bias " + to_string(i) + " doesn't fit word size");
                }
                biasMap[i] = (unsigned int) ((int) scaledBias); // looks stupid, but works (up to 32 bits word size!)
            }
            GenericLut *biasLut = new GenericLut(target, "Bias_Lut"+this->myArguments->getId(), biasMap, this->neuronCounterWidth, this->myArguments->getWeightWordSize());
            addSubComponent(biasLut);
            // in port map
            for (unsigned int i = 0; i < this->neuronCounterWidth; i++) {
                this->vhdl << this->declare("Bias_Lut_i" + to_string(i), 1, false) << " <= "
                           << "Neuron_Counter(" << i << ");" << endl;
                inPortMap(biasLut, "i" + to_string(i), "Bias_Lut_i" + to_string(i));
            }
            // out port map
            for (unsigned int i = 0; i < this->myArguments->getWeightWordSize(); i++) {
                outPortMap(biasLut, "o" + to_string(i), "Bias_" + to_string(i) + "_std", true);
                this->vhdl << "Bias_temp(" << i << ") <= Bias_" << i << "_std;" << endl;
            }
            this->vhdl << instance(biasLut, "Bias_Lut_instance");
        }
        resizeBias();
    }

    void FullConnectedLayer::resizeBias()
    {
        // sign extension
        this->vhdl << declare("Bias_temp_MSB",(this->intermediateResultWidth-this->myArguments->getWeightWordSize()-this->myArguments->getFraction())) << " <= (others => Bias_temp(Bias_temp'length-1));" << endl;
        this->vhdl << "Bias <= Bias_temp_MSB & Bias_temp";

        // fill lower bits
        int loopMax = this->myArguments->getFraction();
        for(int i=0; i<loopMax; i++)
        {
            if(i==0) {this->vhdl << " & \"";}
            this->vhdl << "0";
            if(i==loopMax-1) {this->vhdl << "\"";}
        }
        this->vhdl << ";" << endl;
    }

    void FullConnectedLayer::buildFullConnectedCore(Target* target)
    {
        // full connected core
        FullConnectedCore* core = new FullConnectedCore(target,this->parent->getMultiplicationType(),this->myArguments->getWordSize(),this->myArguments->getFraction(),this->myArguments->getWeightWordSize(),this->myArguments->getWeightFraction(),this->intermediateResultWidth,roundingTypeEnum::truncation);
        addSubComponent(core);
        inPortMap(core,"X",this->getInputSignalName(0));
        inPortMap(core,"Y","Full_Connected_Core_Second_Input");
        inPortMap(core,"W","Full_Connected_Core_Weight_Input");
        outPortMap(core,"R","Full_Connected_Core_Output",false);
        this->vhdl << instance(core, "Full_Connected_Core");

        this->fullConnectedCorePipelineDepth = core->getPipelineDepth();

        // register to delay Intermediate_Result_Reg_enable for the pipeline depth of FullConnectedCore
        Register* regValid = new Register(target,1,this->fullConnectedCorePipelineDepth,false);
        addSubComponent(regValid);
        inPortMap(regValid,"X","Intermediate_Result_Reg_enable");
        outPortMap(regValid,"R","Intermediate_Result_Reg_enable_reg",false);
        this->vhdl << instance(regValid,"ValidData_i_Reg");

        if(this->fullConnectedCorePipelineDepth==0)
        {
            this->vhdl << "Intermediate_Result_Reg_Input <= Full_Connected_Core_Output;" << endl;
        }
        else
        {
            declare("Multiplexer_Switch_temp",1);
            declare("Multiplexer_Switch",1);
            this->buildFullConnectedCoreAdder(target);
        }

        // register to save intermediate result
        Register* reg = new Register(target,this->intermediateResultWidth,1,true,true,false);
        addSubComponent(reg);
        inPortMap(reg,"X","Intermediate_Result_Reg_Input");
        inPortMap(reg,"E","Intermediate_Result_Reg_enable_reg");
        inPortMap(reg,"ExternalReset","Intermediate_Result_Reg_reset");
        outPortMap(reg,"R","Intermediate_Result",false);
        this->vhdl << instance(reg, "Intermediate_Result_Reg");

        // mux to select second input of FullConnectedCore
        this->vhdl << declare("Weight_Counter_Equals_Zero",1) << " <= \"1\" when unsigned(Weight_Counter)=0 else \"0\";" << endl;

        GenericMux* mux = new GenericMux(target,this->intermediateResultWidth,2);
        addSubComponent(mux);
        inPortMap(mux,mux->getSelectName(),"Weight_Counter_Equals_Zero");
        inPortMap(mux,mux->getInputName(0),"Intermediate_Result");
        inPortMap(mux,mux->getInputName(1),"Bias");
        outPortMap(mux,mux->getOutputName(),"Full_Connected_Core_Second_Input",false);
        this->vhdl << instance(mux,"Second_Input_Mux");
    }

    void FullConnectedLayer::buildFullConnectedCoreAdder(Target* target)
    {
        // build registers as adder inputs
        this->vhdl << declare("Full_Connected_Core_Delay_0",this->intermediateResultWidth) << " <= Full_Connected_Core_Output;" << endl;
        if((this->weightsPerNeuron % 32) <= this->fullConnectedCorePipelineDepth)
        {
            declare("Adder_Inputs_External_Reset",1);
            declare("Adder_Inputs_External_Reset_temp",1);
        }

        for(unsigned int i=1; i<=this->fullConnectedCorePipelineDepth; i++)
        {
            declare("Full_Connected_Core_Delay_"+to_string(i),this->intermediateResultWidth);
            Register* reg = new Register(target,this->intermediateResultWidth,1,false,true,true);
            addSubComponent(reg);
            inPortMap(reg,"X","Full_Connected_Core_Delay_"+to_string(i-1));
            if((this->weightsPerNeuron % 32) <= this->fullConnectedCorePipelineDepth)
            {
                inPortMap(reg,"ExternalReset","Adder_Inputs_External_Reset");
            }
            else
            {
                inPortMap(reg,"ExternalReset","Multiplexer_Switch");
            }
            outPortMap(reg,"R","Full_Connected_Core_Delay_"+to_string(i),false);
            this->vhdl << instance(reg,"Full_Connected_Core_Delay_Register_"+to_string(i-1));
        }

        // build adder
        int numberOfInputs = this->fullConnectedCorePipelineDepth+1;
        vector<int> bitheapWIns(numberOfInputs,this->intermediateResultWidth);
        vector<bool> bitheapSigns(numberOfInputs,true);
        vector<int> bitheapWeights(numberOfInputs,0);
        BitheapWrapper* bitH = new BitheapWrapper(target,bitheapWIns,bitheapWIns.size(),bitheapSigns,bitheapWeights);
        addSubComponent(bitH);
        for(unsigned int i=0; i<=this->fullConnectedCorePipelineDepth; i++)
        {
            inPortMap(bitH,"X"+to_string(i),"Full_Connected_Core_Delay_"+to_string(i));
        }
        declare("Full_Connected_Core_Adder_Output",bitH->getOutputWordSize());
        outPortMap(bitH,"R","Full_Connected_Core_Adder_Output",false);
        this->vhdl << instance(bitH,"Full_Connected_Core_Adder_instance");
        this->fullConnectedCoreAdderPipelineDepth = bitH->getPipelineDepth();

        // truncate output of adder
        declare("Full_Connected_Core_Adder_Output_Rounded",this->intermediateResultWidth);
        Rounding* round = new Rounding(target,bitH->getOutputWordSize(),0,this->intermediateResultWidth,0,roundingTypeEnum::truncation);
        addSubComponent(round);
        inPortMap(round,"X","Full_Connected_Core_Adder_Output");
        outPortMap(round,"R","Full_Connected_Core_Adder_Output_Rounded",false);
        this->vhdl << instance(round,"Full_Connected_Core_Adder_Round_instance");

        // build register for Multiplexer_Switch-signal
        Register* reg2 = new Register(target,1,this->fullConnectedCorePipelineDepth,false);
        addSubComponent(reg2);
        inPortMap(reg2,"X","Multiplexer_Switch_temp");
        outPortMap(reg2,"R","Multiplexer_Switch",false);
        this->vhdl << instance(reg2,"Multiplexer_Switch_register_instance");

        // build multiplexer as register input
        GenericMux* multiplexer = new GenericMux(target,this->intermediateResultWidth,2);
        addSubComponent(multiplexer);
        inPortMap(multiplexer,multiplexer->getSelectName(),"Multiplexer_Switch");
        inPortMap(multiplexer,multiplexer->getInputName(0),"Full_Connected_Core_Output");
        inPortMap(multiplexer,multiplexer->getInputName(1),"Full_Connected_Core_Adder_Output_Rounded");
        outPortMap(multiplexer,multiplexer->getOutputName(),"Intermediate_Result_Reg_Input",false);
        this->vhdl << instance(multiplexer,"Multiplexer_For_Register_Input");
    }

    void FullConnectedLayer::buildController(Target* target)
    {
        // controller
        FullConnectedLayerController* contr = new FullConnectedLayerController(target,this->myArguments->getNumberOfOutputFeatures(),this->weightsPerNeuron,this->weightsPerDMAAccess,this->fullConnectedCorePipelineDepth,this->fullConnectedCoreAdderPipelineDepth);
        addSubComponent(contr);
        inPortMap(contr,"newStep","newStep");
        inPortMap(contr,"validData_i",Layer::getValidDataName(0)+"_i");
        inPortMap(contr,"Weights_Are_Valid","Weights_Are_Valid");

        if(contr->thisHasMultiplexerSwitchOutput()==true)
        {
            outPortMap(contr,"Multiplexer_Switch","Multiplexer_Switch_temp",false);
        }
        if((this->weightsPerNeuron % 32) <= this->fullConnectedCorePipelineDepth)
        {
            outPortMap(contr,"Reset_Adder_Inputs","Adder_Inputs_External_Reset_temp",false);
            Register* adderInputsResetReg = new Register(target,1,this->fullConnectedCorePipelineDepth-2);
            addSubComponent(adderInputsResetReg);
            inPortMap(adderInputsResetReg,"X","Adder_Inputs_External_Reset_temp");
            outPortMap(adderInputsResetReg,"R","Adder_Inputs_External_Reset",false);
            vhdl << instance(adderInputsResetReg,"Adder_Inputs_Reset_Reg_instance");
        }

        outPortMap(contr,"finished","finished_temp",false);
        outPortMap(contr,"getNewData",Layer::getGetNewDataName(0),false);
        outPortMap(contr,"getNewData_reset",Layer::getGetNewDataName(0)+"_reset",false);
        outPortMap(contr,"validData_o",Layer::getValidDataName(0)+"_o",false);
        outPortMap(contr,"weightRamEnable","Weight_Register_Enable",false);
        outPortMap(contr,"intermediateResultRegisterEnable","Intermediate_Result_Reg_enable",false);
        outPortMap(contr,"intermediateResultRegisterReset","Intermediate_Result_Reg_reset",false);
        outPortMap(contr,"Weight_Counter","Weight_Counter",false);
        outPortMap(contr,"Neuron_Counter","Neuron_Counter",false);
        outPortMap(contr,"Read_New_Weights","Start_Fetching_New_Weights",false);
        vhdl << instance(contr, "Controller_instance");

        // register to delay finished_temp for the pipeline depth of FullConnectedCore plus pipeline depth of adder plus 1
        Register* reg = new Register(target,1,this->fullConnectedCorePipelineDepth+this->fullConnectedCoreAdderPipelineDepth+1,false);
        addSubComponent(reg);
        inPortMap(reg,"X","finished_temp");
        outPortMap(reg,"R","finished",false);
        this->vhdl << instance(reg,"Finished_Reg");
    }

    void FullConnectedLayer::buildWeightFetcher(Target* target)
    {
        unsigned int memoryAccessesPerNeuron = (unsigned int)ceil((double)this->weightsPerNeuron / (double)this->weightsPerDMAAccess);
        WeightFetcher* fetch = new WeightFetcher(target,this->weightsPerDMAAccess,this->myArguments->getWeightWordSize(),this->myArguments->getStartAddress(),memoryAccessesPerNeuron-1,this->myArguments->getNumberOfOutputFeatures()-1,this->weightsPerNeuron);
        addSubComponent(fetch);

        if(this->weightsPerNeuron>32)
        {
            vhdl << declare("Weight_Counter_MSBs",this->weightCounterWidth-log2(this->weightsPerDMAAccess)) << " <= Weight_Counter(" << this->weightCounterWidth-1 << " downto " << log2(this->weightsPerDMAAccess) << ");" << endl;
            vhdl << declare("Weight_Counter_LSBs",log2(this->weightsPerDMAAccess)) << " <= Weight_Counter(" << log2(this->weightsPerDMAAccess)-1 << " downto " << 0 << ");" << endl;
        }
        else
        {
            vhdl << declare("Weight_Counter_LSBs",log2(this->weightsPerDMAAccess)) << " <= ";
            // fill MSBs of LSB-vector with zeros
            for(unsigned int i=0; i<log2(this->weightsPerDMAAccess)-this->weightCounterWidth; i++)
            {
                if(i==0) vhdl << "\"";
                vhdl << "0";
                if(i==log2(this->weightsPerDMAAccess)-this->weightCounterWidth-1) vhdl << "\" & ";
            }
            // connect actual weight counter to the remaining LSBs of LSB-vector
            vhdl << "Weight_Counter;" << endl;
        }
        inPortMap(fetch,WeightFetcher::getDataInPortName(),Layer::getWeightInputName());
        inPortMap(fetch,WeightFetcher::getValidInPortName(),Layer::getWeightsValidInputName());
        inPortMap(fetch,WeightFetcher::getLastInPortName(),Layer::getLastWeightsInputName());

        if(this->weightsPerNeuron>32)
            inPortMap(fetch,WeightFetcher::getInnerCounterPortName(),"Weight_Counter_MSBs");
        if(this->myArguments->getNumberOfOutputFeatures()>1)
            inPortMap(fetch,WeightFetcher::getOuterCounterPortName(),"Neuron_Counter");

        inPortMap(fetch,WeightFetcher::getStartNewReadPortName(),"Start_Fetching_New_Weights");
        inPortMap(fetch,WeightFetcher::getWeightsArrivedAtConvCorePortName(),"Weight_Register_Enable");

        outPortMap(fetch,WeightFetcher::getNewReadAccessPortName(),Layer::getNewReadAccessOutputName(),false);
        outPortMap(fetch,WeightFetcher::getNumberOfBytesPortName(),Layer::getNumberOfBytesOutputName(),false);
        outPortMap(fetch,WeightFetcher::getStartAddressPortName(),Layer::getNextStartAddressOutputName(),false);
        outPortMap(fetch,WeightFetcher::getNextWeightsValidPortName(),"Weights_Are_Valid",false);
        for(unsigned int i=0; i<this->weightsPerDMAAccess; i++)
        {
            outPortMap(fetch,fetch->getNextWeightsPortName(i),"Next_Weights_"+to_string(i),true);
        }
        vhdl << instance(fetch,"Weight_Fetcher_instance");
    }

    void FullConnectedLayer::buildWeightRegisters(Target* target)
    {
        // build registers
        for(unsigned int i=0; i<this->weightsPerDMAAccess; i++)
        {
            Register* reg = new Register(target,this->myArguments->getWeightWordSize(),1,true);
            addSubComponent(reg);
            inPortMap(reg,"X","Next_Weights_"+to_string(i));
            inPortMap(reg,"E","Weight_Register_Enable");
            outPortMap(reg,"R","Next_Weights_"+to_string(i)+"_reg",true);
            vhdl << instance(reg, "Weight_Register_instance_"+to_string(i));
        }

        // build mux to select right register
        GenericMux* weightMux = new GenericMux(target,this->myArguments->getWeightWordSize(),this->weightsPerDMAAccess);
        addSubComponent(weightMux);
        inPortMap(weightMux,weightMux->getSelectName(),"Weight_Counter_LSBs");
        for(unsigned int i=0; i<this->weightsPerDMAAccess; i++)
        {
            inPortMap(weightMux,weightMux->getInputName(i),"Next_Weights_"+to_string(i)+"_reg");
        }
        outPortMap(weightMux,weightMux->getOutputName(),"Full_Connected_Core_Weight_Input",false);
        this->vhdl << instance(weightMux,"Weight_Select_Mux_instance");
    }

    void FullConnectedLayer::roundOutput(Target* target)
    {
        flopoco::roundingTypeEnum rType = flopoco::roundingTypeEnum::saturation;
        if(this->roundingType==0x00) rType = flopoco::roundingTypeEnum::truncation;
        Rounding* round = new Rounding(target,this->intermediateResultWidth,this->myArguments->getFraction()+this->myArguments->getWeightFraction(),this->myArguments->getWordSize(),this->myArguments->getFraction(),rType);
        addSubComponent(round);
        inPortMap(round,"X","Intermediate_Result");
        outPortMap(round,"R","R_rounded",false);
        vhdl << instance(round,"Rounding_instance");
    }

    void FullConnectedLayer::buildActivationFunction(Target* target)
    {
        if(this->myArguments->getActivationFunction()=="relu")
        {
            ReLU* relu = new ReLU(target,this->myArguments->getWordSize());
            addSubComponent(relu);
            inPortMap(relu,"X","R_rounded");
            outPortMap(relu,"R",this->getOutputSignalName(0),false);
            vhdl << instance(relu,"ReLU_instance");
        }
        else
        {
            vhdl << this->getOutputSignalName(0) << " <= R_rounded;" << endl;
        }
    }


}//namespace flopoco
