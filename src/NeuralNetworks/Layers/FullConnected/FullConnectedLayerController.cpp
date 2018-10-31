// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "FullConnectedLayerController.hpp"

#include "NeuralNetworks/Utility/BitheapWrapper.hpp"
#include "NeuralNetworks/Utility/Rounding.hpp"
#include "NeuralNetworks/Utility/Register.hpp"

using namespace std;
namespace flopoco {

    FullConnectedLayerController::FullConnectedLayerController(Target* target, unsigned int numberOfNeurons, unsigned int weightsPerNeuron, unsigned int weightRAMDepth, unsigned int fullConnectedCorePipelineDepth, int fullConnectedCoreAdderPipelineDepth)
            : Operator(target)
    {
        useNumericStd();
        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="FullConnectedLayerController";

        // author
        setCopyrightString("Nicolai Fiege, 2018");

        // definition of the name of the operator
        ostringstream name;
        name << "FullConnectedLayerController_" << numberOfNeurons << "_" << weightsPerNeuron << "_" << weightRAMDepth << "_" << fullConnectedCorePipelineDepth << "_" << ((fullConnectedCoreAdderPipelineDepth<0)?("minus1"):(to_string(fullConnectedCoreAdderPipelineDepth)));
        setName(name.str());

        unsigned int neuronCounterWidth = ceil(log2(numberOfNeurons));
        unsigned int weightCounterWidth = ceil(log2(weightsPerNeuron));

        // weightRAMDepth must be a power of 2
        unsigned int weightRAMAddressWidth = ceil(log2(weightRAMDepth));
        if(weightRAMDepth != ((weightRAMDepth >> weightRAMAddressWidth) << weightRAMAddressWidth) || weightRAMAddressWidth == 0)
        {
            THROWERROR("Weight RAM Depth must be a power of 2 and bigger than 1");
        }

        // write vhdl code

        ///////////
        // ports //
        ///////////
        addInput("newStep",1);
        addOutput("finished",1);

        addInput("validData_i",1);
        addOutput("getNewData",1);
        addOutput("getNewData_reset",1);
        addOutput("validData_o",1);

        addOutput("weightRamEnable",1);

        addOutput("intermediateResultRegisterEnable",1);
        addOutput("intermediateResultRegisterReset",1);

        addOutput("Weight_Counter",weightCounterWidth);
        addOutput("Neuron_Counter",neuronCounterWidth);

        addOutput("Read_New_Weights",1);
        addInput("Weights_Are_Valid",1);

        /////////////////////
        // declare signals //
        /////////////////////
        declare("validData_o_temp",1);
        declare("validData_o_temp_reg",1);
        declare("validData_o_temp_reg2",1);
        declare("Read_New_Weights_temp",1);

        this->hasMultiplexerSwitchOutput = (fullConnectedCoreAdderPipelineDepth >= 0);

        if(fullConnectedCoreAdderPipelineDepth<0)
        {
            /////////////////////////////////////////
            // enable intermediate result register //
            /////////////////////////////////////////
            this->vhdl << "intermediateResultRegisterEnable <= validData_i;" << endl;
        }
        else
        {
            /////////////////////////////////////////
            // enable intermediate result register //
            /////////////////////////////////////////
            this->vhdl << "intermediateResultRegisterEnable <= validData_i or " << declare("regEnableFromCounterValue",1) << ";" << endl;
            if(weightCounterWidth >= weightRAMAddressWidth)
            {
                this->vhdl << declare("regEnableFromCounterValue_temp",1) << " <= \"1\" when (unsigned(Weight_Counter_temp("
                           << log2(weightRAMDepth)-1 << " downto 0)) = " << weightRAMDepth-1
                           << " or unsigned(Weight_Counter_temp)=" << weightsPerNeuron-1 << ") else \"0\";" << endl;
            }
            else
            {
                this->vhdl << declare("regEnableFromCounterValue_temp",1) << " <= \"1\" when unsigned(Weight_Counter_temp)="
                           << weightsPerNeuron-1 << " else \"0\";" << endl;
            }

            Register* regEnableReg = new Register(target,1,fullConnectedCoreAdderPipelineDepth,false);
            addSubComponent(regEnableReg);
            inPortMap(regEnableReg,"X","regEnableFromCounterValue_temp");
            outPortMap(regEnableReg,"R","regEnableFromCounterValue",false);
            this->vhdl << instance(regEnableReg,"Reg_Enable_From_Counter_Value_Register_instance");

            ////////////////////////////////////////////////////////////////////////////////////////////
            // Generate switch-signal for Mux (that selects the input to the intermediate result reg) //
            ////////////////////////////////////////////////////////////////////////////////////////////
            addOutput("Multiplexer_Switch",1);
            this->vhdl << "Multiplexer_Switch <= regEnableFromCounterValue;" << endl;

            ////////////////////////////////////////
            // reset intermediate result register //
            ////////////////////////////////////////
            if(fullConnectedCorePipelineDepth<2) THROWERROR("Need full connected core to have pipeline depth > 1");
            if(weightCounterWidth >= weightRAMAddressWidth)
            {
                this->vhdl << declare("intermediateResultRegResetFromCounterValues",1)
                           << " <= \"1\" when unsigned(Weight_Counter_temp(" << weightRAMAddressWidth-1
                           << " downto 0)) >= 1 and unsigned(Weight_Counter_temp(" << weightRAMAddressWidth-1
                           << " downto 0)) <= " << fullConnectedCorePipelineDepth-1 << " else \"0\";" << endl;
            }
            else
            {
                this->vhdl << declare("intermediateResultRegResetFromCounterValues",1)
                           << " <= \"1\" when unsigned(Weight_Counter_temp) >= 1 and unsigned(Weight_Counter_temp) <= "
                           << fullConnectedCorePipelineDepth-1 << " else \"0\";" << endl;
            }

            Register* oneMoreDelay = new Register(target,1,1,false,false); // need additional delay, because the register has an async reset
            addSubComponent(oneMoreDelay);
            inPortMap(oneMoreDelay,"X","validData_o_temp_reg");
            outPortMap(oneMoreDelay,"R","validData_o_temp_reg2",false);
            this->vhdl << instance(oneMoreDelay,"ValidData_o_Reg_Delay1");
            this->vhdl << "intermediateResultRegisterReset <= validData_o_temp_reg2 or intermediateResultRegResetFromCounterValues;" << endl;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // reset adder inputs if the number of weights in the last package is smaller than the pipeline depth of the fullConnectedCore //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(weightsPerNeuron % 32 <= fullConnectedCorePipelineDepth)
        {
            declare("validData_i_reg",1);
            addOutput("Reset_Adder_Inputs",1);
            Register* validReg = new Register(target,1,1);
            addSubComponent(validReg);
            inPortMap(validReg,"X","validData_i");
            outPortMap(validReg,"R","validData_i_reg",false);
            this->vhdl << instance(validReg,"ValidData_i_Reg_instance");

            this->vhdl << "process(clk)" << endl;
            this->vhdl << "begin" << endl;
            this->vhdl << tab << "if(rising_edge(clk)) then" << endl;
            this->vhdl << tab << tab << "if(rst='1') then" << endl;
            this->vhdl << tab << tab << tab << "Reset_Adder_Inputs <= \"0\";" << endl;
            this->vhdl << tab << tab << "elsif(validData_i=\"1\" and validData_i_reg=\"0\") then" << endl;
            this->vhdl << tab << tab << tab << "Reset_Adder_Inputs <= \"1\";" << endl;
            this->vhdl << tab << tab << "else" << endl;
            this->vhdl << tab << tab << tab << "Reset_Adder_Inputs <= \"0\";" << endl;
            this->vhdl << tab << tab << "end if;" << endl;
            this->vhdl << tab << "end if;" << endl;
            this->vhdl << "end process;" << endl;
        }

        //////////////
        // counters //
        //////////////
        this->vhdl << "Weight_Counter <= " << declare("Weight_Counter_temp",weightCounterWidth) << ";" << endl;
        this->vhdl << "Neuron_Counter <= " << declare("Neuron_Counter_temp",neuronCounterWidth) << ";" << endl;
        this->vhdl << "process(clk)" << endl;
        this->vhdl << "begin" << endl;
        this->vhdl << tab << "if(rising_edge(clk)) then" << endl;
        this->vhdl << tab << tab << "if(rst='1' or newStep=\"1\") then" << endl;
        this->vhdl << tab << tab << tab << "Weight_Counter_temp <= (others => '0');" << endl;
        this->vhdl << tab << tab << tab << "Neuron_Counter_temp <= (others => '0');" << endl;
        this->vhdl << tab << tab << "elsif(validData_i=\"1\") then" << endl;
        this->vhdl << tab << tab << tab << "if(unsigned(Weight_Counter_temp) = " << weightsPerNeuron-1 << ") then" << endl;
        this->vhdl << tab << tab << tab << tab << "Weight_Counter_temp <= (others => '0');" << endl;
        this->vhdl << tab << tab << tab << tab << "if(unsigned(Neuron_Counter_temp) = " << numberOfNeurons-1 << ") then" << endl;
        this->vhdl << tab << tab << tab << tab << tab << "Neuron_Counter_temp <= (others => '0');" << endl;
        this->vhdl << tab << tab << tab << tab << "else" << endl;
        this->vhdl << tab << tab << tab << tab << tab << "Neuron_Counter_temp <= std_logic_vector(unsigned(Neuron_Counter_temp) + 1);" << endl;
        this->vhdl << tab << tab << tab << tab << "end if;" << endl;
        this->vhdl << tab << tab << tab << "else" << endl;
        this->vhdl << tab << tab << tab << tab << "Weight_Counter_temp <= std_logic_vector(unsigned(Weight_Counter_temp) + 1);" << endl;
        this->vhdl << tab << tab << tab << "end if;" << endl;
        this->vhdl << tab << tab << "end if;" << endl;
        this->vhdl << tab << "end if;" << endl;
        this->vhdl << "end process;" << endl;
		
		//////////////////////////////////////////////////////////////////////
		// this layer finished calculation and waits for the newStep signal //
		//////////////////////////////////////////////////////////////////////
		this->vhdl << "finished <= " << declare("finished_temp",1) << ";" << endl;
		this->vhdl << "process(clk)" << endl;
		this->vhdl << "begin" << endl;
		this->vhdl << tab << "if(rising_edge(clk)) then" << endl;
		this->vhdl << tab << tab << "if(rst = '1' or newStep = \"1\") then" << endl;
		this->vhdl << tab << tab << tab << "finished_temp <= \"0\";" << endl;
		this->vhdl << tab << tab << "elsif(validData_i=\"1\" and unsigned(Weight_Counter_temp) = " << weightsPerNeuron-1 << " and unsigned(Neuron_Counter_temp) = " << numberOfNeurons-1 << ") then" << endl;
		this->vhdl << tab << tab << tab << "finished_temp <= \"1\";" << endl;
		this->vhdl << tab << tab << "end if;" << endl;
		this->vhdl << tab << "end if;" << endl;
		this->vhdl << "end process;" << endl;

        /////////////////////////////////////////
        // manage reading new weights with DMA //
        /////////////////////////////////////////
        declare("Weights_Arrived",1);

        vhdl << declare("Need_New_Weights_From_Weights_Counter_Value",1) << "(0) <= (Weight_Counter_temp(0)";
        if(weightCounterWidth >= weightRAMAddressWidth)
        {
            for(unsigned int i=1; i<weightRAMAddressWidth; i++)
            {
                vhdl << " and Weight_Counter_temp(" << i << ")";
            }
        }
        else
        {
            for(unsigned int i=1; i<weightCounterWidth; i++)
            {
                vhdl << " and Weight_Counter_temp(" << i << ")";
            }
        }
        vhdl << ") or (not Weights_Counter_Smaller_Max(0));" << endl;


        vhdl << "process(clk)" << endl;
        vhdl << "begin" << endl;
        vhdl << tab << "if(rising_edge(clk)) then" << endl;
        vhdl << tab << tab << "if(rst = '1' or newStep = \"1\" or Need_New_Weights_From_Weights_Counter_Value = \"1\") then" << endl;
        vhdl << tab << tab << tab << "Weights_Arrived <= \"0\";" << endl;
        vhdl << tab << tab << tab << "Read_New_Weights_temp <= \"1\";" << endl;
        vhdl << tab << tab << "elsif(Weights_Arrived=\"0\" and Weights_Are_Valid=\"1\") then" << endl;
        vhdl << tab << tab << tab << "Weights_Arrived <= \"1\";" << endl;
        vhdl << tab << tab << tab << "Read_New_Weights_temp <= \"0\";" << endl;
        vhdl << tab << tab << "end if;" << endl;
        vhdl << tab << "end if;" << endl;
        vhdl << "end process;" << endl;

        vhdl << "process(clk)" << endl;
        vhdl << "begin" << endl;
        vhdl << tab << "if(rising_edge(clk)) then" << endl;
        vhdl << tab << tab << "if(rst = '1' or newStep = \"1\" or Need_New_Weights_From_Weights_Counter_Value = \"1\" or finished_temp=\"1\") then" << endl;
        vhdl << tab << tab << tab << "Read_New_Weights <= \"0\";" << endl;
        vhdl << tab << tab << "else" << endl;
        vhdl << tab << tab << tab << "Read_New_Weights <= Read_New_Weights_temp;" << endl;
        vhdl << tab << tab << "end if;" << endl;
        vhdl << tab << "end if;" << endl;
        vhdl << "end process;" << endl;

        this->vhdl << "weightRamEnable <= not Weights_Arrived and Weights_Are_Valid;" << endl;

        ///////////////////////////////////////
        // manage getting new data from BRAM //
        ///////////////////////////////////////
        // => read new weights as long as the current weights are valid and the counters have the "right" values
		vhdl << declare("Weights_Counter_Smaller_Max",1) << " <= \"1\" when unsigned(Weight_Counter_temp) < " << weightsPerNeuron-1 << " else \"0\";" << endl;
		vhdl << declare("Neuron_Counter_Smaller_Max",1) << " <= \"1\" when unsigned(Neuron_Counter_temp) < " << numberOfNeurons-1 << " else \"0\";" << endl;
        vhdl << "getNewData <= (not finished_temp) and Weights_Arrived and (not Need_New_Weights_From_Weights_Counter_Value) and (Neuron_Counter_Smaller_Max or Weights_Counter_Smaller_Max);" << endl;

        vhdl << "process(clk)" << endl;
        vhdl << "begin" << endl;
        vhdl << tab << "if(rising_edge(clk)) then" << endl;
        vhdl << tab << tab << "if(rst='1') then" << endl;
        vhdl << tab << tab << tab << "getNewData_reset <= \"0\";" << endl;
        vhdl << tab << tab << "elsif(validData_i=\"1\" and unsigned(Weight_Counter_temp) = " << weightsPerNeuron-1 << ") then" << endl;
        vhdl << tab << tab << tab << "getNewData_reset <= \"1\";" << endl;
        vhdl << tab << tab << "else" << endl;
        vhdl << tab << tab << tab << "getNewData_reset <= \"0\";" << endl;
        vhdl << tab << tab << "end if;" << endl;
        vhdl << tab << "end if;" << endl;
        vhdl << "end process;" << endl;

        ///////////////////////////////////////
        // manage writing valid data to BRAM //
        ///////////////////////////////////////
        vhdl << "validData_o <= validData_o_temp_reg;" << endl;
        vhdl << "process(clk)" << endl;
        vhdl << "begin" << endl;
        vhdl << tab << "if(rising_edge(clk)) then" << endl;
        vhdl << tab << tab << "if(rst='1' or validData_o_temp = \"1\") then" << endl;
        vhdl << tab << tab << tab << "validData_o_temp <= \"0\";" << endl;
        vhdl << tab << tab << "elsif(validData_i=\"1\" and unsigned(Weight_Counter_temp) = " << weightsPerNeuron-1 << ") then" << endl;
        vhdl << tab << tab << tab << "validData_o_temp <= \"1\";" << endl;
        vhdl << tab << tab << "end if;" << endl;
        vhdl << tab << "end if;" << endl;
        vhdl << "end process;" << endl;

        Register* validDataReg = new Register(target,1,fullConnectedCorePipelineDepth+fullConnectedCoreAdderPipelineDepth+1,false);
        addSubComponent(validDataReg);
        inPortMap(validDataReg,"X","validData_o_temp");
        outPortMap(validDataReg,"R","validData_o_temp_reg",false);
        this->vhdl << instance(validDataReg,"ValidData_o_Reg");
    }
}