// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "OutputLayer.hpp"
#include "FindIndexOfMaximum.hpp"

#include "NeuralNetworks/MemoryManagement/BlockRam.hpp"
#include "NeuralNetworks/Utility/DataGuard.hpp"
#include "NeuralNetworks/Utility/Register.hpp"
#include "NeuralNetworks/Utility/ModuloCounter.hpp"

using namespace std;
namespace flopoco {




    OutputLayer::OutputLayer(Target* target, NeuralNetwork* parent, int numberOfFeatures_, int wordSize_, bool bramOnOutput_, int width_, int height_, bool threeBRAMS_, bool needIntermediateValue_, bool hasSpecificAddressOnResetPort, bool onlyOutputClassIndex) :
        Layer(target,parent), numberOfFeatures(numberOfFeatures_), wordSize(wordSize_), bramOnOutput(bramOnOutput_), width(width_), height(height_), threeBRAMS(threeBRAMS_), needIntermediateValue(needIntermediateValue_) {

        this->myMaxPtr = nullptr;

        if(wordSize<1)
        {
            stringstream e;
            e << "Need word size > 0!";
            THROWERROR(e.str());
        }
        if(numberOfFeatures<1)
        {
            stringstream e;
            e << "Need number of features > 0!";
            THROWERROR(e.str());
        }
        if(bramOnOutput==true)
        {
            if(this->width<1)
            {
                stringstream e;
                e << "Need width > 0 if u want to write data in BRAM!";
                THROWERROR(e.str());
            }
            if(this->height<1)
            {
                stringstream e;
                e << "Need height > 0 if u want to write data in BRAM!";
                THROWERROR(e.str());
            }
        }


        useNumericStd();

		// definition of the source file name, used for info and error reporting using REPORT 
		srcFileName="OutputLayer";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

		// definition of the name of the operator
		ostringstream name;
        name << "OutputLayer";
        setName(name.str());

        addInput("newStep",1);

        if(this->bramOnOutput==false)
        {
            addInput(Layer::getValidDataName(0)+"_i",1);
            addInput("X_"+to_string(0),wordSize);
            // add In/Out
            for(int i=0; i<numberOfFeatures; i++)
            {
                declare("R_"+to_string(i)+"_temp",wordSize);
                if(onlyOutputClassIndex==false)
                {
                    addOutput("R_"+to_string(i),wordSize);
                    vhdl << "R_" << i << " <= R_" << i << "_temp;" << endl;
                }
            }

            // mod counter
            ModuloCounter* modC = new ModuloCounter(target,this->numberOfFeatures,false);
            addSubComponent(modC);
            inPortMap(modC,"enable",Layer::getValidDataName(0)+"_i");
            inPortMap(modC,"manualReset","newStep");
            outPortMap(modC,"counter","Neuron_Counter",true);
            vhdl << instance(modC,"Neuron_Counter_instance");

            for(int i=0; i<numberOfFeatures; i++)
            {
                // DataGuard for Register-enable
                DataGuard* datG = new DataGuard(target,1,getSignalByName("Neuron_Counter")->width(),i);
                addSubComponent(datG);
                inPortMap(datG,"Data_in",Layer::getValidDataName(0)+"_i");
                inPortMap(datG,"Guard_in","Neuron_Counter");
                outPortMap(datG,"Data_out","Register_enable_"+to_string(i),true);
                vhdl << instance(datG,"Register_Enable_Guard_instance_"+to_string(i));

                // register for each output neuron
                Register* reg = new Register(target,this->wordSize,1,true);
                addSubComponent(reg);
                inPortMap(reg,"X","X_"+to_string(0));
                inPortMap(reg,"E","Register_enable_"+to_string(i));
                outPortMap(reg,"R","R_"+to_string(i)+"_temp",false);
                vhdl << instance(reg,"Data_Reg_instance_"+to_string(i));
            }

            // calc index of maximum
            FindIndexOfMaximum* indMax = new FindIndexOfMaximum(target,wordSize,numberOfFeatures);
            addSubComponent(indMax);
            for(int i=0; i<numberOfFeatures; i++)
            {
                inPortMap(indMax,"X"+to_string(i),"R_"+to_string(i)+"_temp");
            }
            outPortMap(indMax,"R","I",true);
            vhdl << instance(indMax,"Max_Index_instance");
            this->myMaxPtr = indMax;

            syncCycleFromSignal("I");

            // enabled register to save the class index when a new step is initiated
            vhdl << declare("Outputs_Are_Valid_temp",1) << " <= newStep;" << endl;

            addOutput("I_max",indMax->getIndexWordSize());
            Register* classIndexReg = new Register(target,indMax->getIndexWordSize(),1,true);
            addSubComponent(classIndexReg);
            inPortMap(classIndexReg,"X","I");
            inPortMap(classIndexReg,"E","Outputs_Are_Valid_temp");
            outPortMap(classIndexReg,"R","I_max",false);
            this->vhdl << instance(classIndexReg,"Class_Index_Register_instance");

            // register to add one additional pipeline stage to the newStep-signal
            addOutput("Outputs_Are_Valid",1);
            Register* newStepReg = new Register(target,1,1);
            addSubComponent(newStepReg);
            inPortMap(newStepReg,"X","Outputs_Are_Valid_temp");
            outPortMap(newStepReg,"R","Outputs_Are_Valid",false);
            this->vhdl << instance(newStepReg,"Outputs_Are_Valid_Register_instance");
        }
        else
        {
            for(int i=0; i<numberOfFeatures; i++)
            {
                addInput("X_"+to_string(i),wordSize);
                addOutput("R_"+to_string(i),wordSize);
                if(needIntermediateValue==true)
                {
                    addOutput(this->getIntermediateResultName(i),wordSize);
                    addInput(this->getGetNewDataName(i)+"_intermediate",1);
                    addInput(this->getGetNewDataName(i)+"_intermediate_reset",1);
                }
                addInput(Layer::getValidDataName(i)+"_i",1);
                addOutput(Layer::getValidDataName(i)+"_o",1);

                if(threeBRAMS==true)
                {
                    // Signals from last layer: we have X_in, validData_in

                    // External Read side: we have R_out
                    addInput("getNewData_"+to_string(i),1);
                    addInput("getNewDataSet_"+to_string(i),1);

                    // flag signals to decide which memory area to read and which one to write
                    // if the last layer is Pooling/Convolutional we need 3 BRAMs per feature
                    // so the writing layer and the external output can access the memory arbitrarily
                    declare("writeFrame_"+to_string(i),2); // can be "00", "01", or "10"
                    declare("readFrame_"+to_string(i),2); // can be "00", "01", or "10"

                    // communication with BRAM
                    int numberOfPixels = this->width * this->height;
                    int BRAMWidth = ceil(log2(numberOfPixels))+2; // +2 for the frame-counters
                    declare("WriteAddress_"+to_string(i),BRAMWidth);
                    declare("ReadAddress_"+to_string(i),BRAMWidth);
                    declare("writeCounter_"+to_string(i),(BRAMWidth-2));
                    declare("readCounter_"+to_string(i),(BRAMWidth-2));
                    if(needIntermediateValue==true)
                    {
                        declare("readIntermediateCounter_"+to_string(i),(BRAMWidth-2));
                        declare("ReadIntermediateAddress_"+to_string(i),BRAMWidth);
                    }
                    if(hasSpecificAddressOnResetPort==true)
                    {
                        addInput("reset_address"+to_string(i),(BRAMWidth-2));
                    }

                    // inst BRAM
                    BlockRam* bram = new BlockRam(target,this->wordSize,BRAMWidth,needIntermediateValue);
                    addSubComponent(bram);
                    inPortMap(bram,"Address_W_in","WriteAddress_"+to_string(i));
                    inPortMap(bram,"Address_R_in","ReadAddress_"+to_string(i));
                    inPortMap(bram,"Data_in","X_"+to_string(i));
                    inPortMap(bram,"WriteEnable",Layer::getValidDataName(i)+"_i");
                    outPortMap(bram,"Data_out","R_"+to_string(i),false);
                    if(needIntermediateValue==true)
                    {
                        outPortMap(bram,"Data2_out",this->getIntermediateResultName(i),false);
                        inPortMap(bram,"Address_R2_in","ReadIntermediateAddress_"+to_string(i));
                    }

                    this->vhdl << instance(bram,"BRAM_"+to_string(i));

                    this->vhdl << "WriteAddress_" << i << " <= writeFrame_" << i << " & writeCounter_" << i << ";" << endl;
                    this->vhdl << "ReadAddress_" << i << " <= readFrame_" << i << " & readCounter_" << i << ";" << endl << endl;
                    if(needIntermediateValue==true)
                    {
                        this->vhdl << "ReadIntermediateAddress_" << i << " <= writeFrame_" << i << " & readIntermediateCounter_" << i << ";" << endl << endl;
                    }

                    // counters
                    this->vhdl << "process(clk)" << endl;
                    this->vhdl << "begin" << endl;
                    this->vhdl << tab << "if(rising_edge(clk)) then" << endl;
                    this->vhdl << tab << tab << "if(rst='1') then" << endl;
                    this->vhdl << tab << tab << tab << "writeCounter_" << i << " <= (others => '0');" << endl;
                    this->vhdl << tab << tab << tab << "readCounter_" << i << " <= (others => '0');" << endl;
                    if(needIntermediateValue==true)
                    {
                        this->vhdl << tab << tab << tab << "readIntermediateCounter_" << i << " <= (others => '0');" << endl;
                    }
                    this->vhdl << tab << tab << tab << Layer::getValidDataName(i)+"_o" << " <= \"0\";" << endl;
                    this->vhdl << tab << tab << "else" << endl;
                    if(hasSpecificAddressOnResetPort==true)
                    {
                        this->vhdl << tab << tab << tab << "if(newStep=\"1\") then" << endl;
                        this->vhdl << tab << tab << tab << tab << "writeCounter_" << i << " <= (others => '0');" << endl;
                        this->vhdl << tab << tab << tab << tab << "readIntermediateCounter_" << i << " <= (others => '0');" << endl;
                        this->vhdl << tab << tab << tab << "elsif(" << this->getGetNewDataName(i) << "_intermediate_reset=\"1\") then" << endl;
                        this->vhdl << tab << tab << tab << tab << "writeCounter_" << i << " <= reset_address;" << endl;
                        this->vhdl << tab << tab << tab << tab << "readIntermediateCounter_" << i << " <= reset_address;" << endl;
                    }
                    else
                    {
                        this->vhdl << tab << tab << tab << "if(newStep=\"1\"" << ((needIntermediateValue==true)?(" or "+this->getGetNewDataName(i)+"_intermediate_reset"+" = \"1\""):("")) << ") then" << endl;
                        this->vhdl << tab << tab << tab << tab << "writeCounter_" << i << " <= (others => '0');" << endl;
                        if(needIntermediateValue==true)
                        {
                            this->vhdl << tab << tab << tab << tab << "readIntermediateCounter_" << i << " <= (others => '0');" << endl;
                        }
                    }
                    this->vhdl << tab << tab << tab << "else" << endl;
                    this->vhdl << tab << tab << tab << tab << "if(" << Layer::getValidDataName(i)+"_i" << "=\"1\") then" << endl;
                    this->vhdl << tab << tab << tab << tab << tab << "writeCounter_" << i << " <= std_logic_vector(unsigned(writeCounter_" << i << ")+1);" << endl;
                    this->vhdl << tab << tab << tab << tab << "end if;" << endl;
                    if(needIntermediateValue==true)
                    {
                        this->vhdl << tab << tab << tab << tab << "if(" << this->getGetNewDataName(i) + "_intermediate" << "=\"1\") then" << endl;
                        this->vhdl << tab << tab << tab << tab << tab << "readIntermediateCounter_" << i << " <= std_logic_vector(unsigned(readIntermediateCounter_" << i << ")+1);" << endl;
                        this->vhdl << tab << tab << tab << tab << "end if;" << endl;
                    }
                    this->vhdl << tab << tab << tab << "end if;" << endl;
                    this->vhdl << tab << tab << tab << "if(getNewDataSet_" << i << "=\"1\") then" << endl;
                    this->vhdl << tab << tab << tab << tab << "readCounter_" << i << " <= (others => '0');" << endl;
                    this->vhdl << tab << tab << tab << tab << Layer::getValidDataName(i)+"_o" << " <= \"0\";" << endl;
                    this->vhdl << tab << tab << tab << "else" << endl;
                    this->vhdl << tab << tab << tab << tab << "if(getNewData_" << i << "=\"1\") then" << endl;
                    this->vhdl << tab << tab << tab << tab << tab << "readCounter_" << i << " <= std_logic_vector(unsigned(readCounter_" << i << ")+1);" << endl;
                    this->vhdl << tab << tab << tab << tab << tab << Layer::getValidDataName(i)+"_o" << " <= \"1\";" << endl;
                    this->vhdl << tab << tab << tab << tab << "else" << endl;
                    this->vhdl << tab << tab << tab << tab << tab << Layer::getValidDataName(i)+"_o" << " <= \"0\";" << endl;
                    this->vhdl << tab << tab << tab << tab << "end if;" << endl;
                    this->vhdl << tab << tab << tab << "end if;" << endl;
                    this->vhdl << tab << tab << "end if;" << endl;
                    this->vhdl << tab << "end if;" << endl;
                    this->vhdl << "end process;" << endl << endl;

                    // frame counters
                    this->vhdl << "process(clk)" << endl
                               << "begin" << endl
                               << tab << "if(rising_edge(clk)) then" << endl
                               << tab << tab << "if(rst='1') then" << endl
                               << tab << tab << tab << "writeFrame_" << i << " <= \"00\";" << endl
                               << tab << tab << tab << "readFrame_" << i << " <= \"01\";" << endl
                               << tab << tab << "else" << endl
                               << tab << tab << tab << "if(newStep=\"1\") then" << endl
                               << tab << tab << tab << tab << "case writeFrame_" << i << " is" << endl
                               << tab << tab << tab << tab << tab << "when \"00\" =>" << endl
                               << tab << tab << tab << tab << tab << tab << "case readFrame_" << i << " is" << endl
                               << tab << tab << tab << tab << tab << tab << tab << "when \"01\" => writeFrame_" << i << " <= \"10\";" << endl
                               << tab << tab << tab << tab << tab << tab << tab << "when others => writeFrame_" << i << " <= \"01\";" << endl
                               << tab << tab << tab << tab << tab << tab << "end case;" << endl
                               << tab << tab << tab << tab << tab << "when \"01\" =>" << endl
                               << tab << tab << tab << tab << tab << tab << "case readFrame_" << i << " is" << endl
                               << tab << tab << tab << tab << tab << tab << tab << "when \"00\" => writeFrame_" << i << " <= \"10\";" << endl
                               << tab << tab << tab << tab << tab << tab << tab << "when others => writeFrame_" << i << " <= \"00\";" << endl
                               << tab << tab << tab << tab << tab << tab << "end case;" << endl
                               << tab << tab << tab << tab << tab << "when others =>" << endl
                               << tab << tab << tab << tab << tab << tab << "case readFrame_" << i << " is" << endl
                               << tab << tab << tab << tab << tab << tab << tab << "when \"00\" => writeFrame_" << i << " <= \"01\";" << endl
                               << tab << tab << tab << tab << tab << tab << tab << "when others => writeFrame_" << i << " <= \"00\";" << endl
                               << tab << tab << tab << tab << tab << tab << "end case;" << endl
                               << tab << tab << tab << tab << "end case;" << endl
                               << tab << tab << tab << "end if;" << endl << endl
                               << tab << tab << tab << "if(getNewDataSet_" << i << "=\"1\") then" << endl
                               << tab << tab << tab << tab << "case readFrame_" << i << " is" << endl
                               << tab << tab << tab << tab << tab << "when \"00\" =>" << endl
                               << tab << tab << tab << tab << tab << tab << "case writeFrame_" << i << " is" << endl
                               << tab << tab << tab << tab << tab << tab << tab << "when \"01\" =>" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "if(newStep=\"1\") then" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << tab << "readFrame_" << i << " <= \"01\";" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "else" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << tab << "readFrame_" << i << " <= \"10\";" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "end if;" << endl
                               << tab << tab << tab << tab << tab << tab << tab << "when others => " << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "if(newStep=\"1\") then" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << tab << "readFrame_" << i << " <= \"10\";" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "else" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << tab << "readFrame_" << i << " <= \"01\";" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "end if;" << endl
                               << tab << tab << tab << tab << tab << tab << "end case;" << endl
                               << tab << tab << tab << tab << tab << "when \"01\" =>" << endl
                               << tab << tab << tab << tab << tab << tab << "case writeFrame_" << i << " is" << endl
                               << tab << tab << tab << tab << tab << tab << tab << "when \"00\" =>" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "if(newStep=\"1\") then" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << tab << "readFrame_" << i << " <= \"00\";" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "else" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << tab << "readFrame_" << i << " <= \"10\";" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "end if;" << endl
                               << tab << tab << tab << tab << tab << tab << tab << "when others => " << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "if(newStep=\"1\") then" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << tab << "readFrame_" << i << " <= \"10\";" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "else" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << tab << "readFrame_" << i << " <= \"00\";" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "end if;" << endl
                               << tab << tab << tab << tab << tab << tab << "end case;" << endl
                               << tab << tab << tab << tab << tab << "when others =>" << endl
                               << tab << tab << tab << tab << tab << tab << "case writeFrame_" << i << " is" << endl
                               << tab << tab << tab << tab << tab << tab << tab << "when \"00\" =>" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "if(newStep=\"1\") then" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << tab << "readFrame_" << i << " <= \"00\";" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "else" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << tab << "readFrame_" << i << " <= \"01\";" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "end if;" << endl
                               << tab << tab << tab << tab << tab << tab << tab << "when others => " << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "if(newStep=\"1\") then" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << tab << "readFrame_" << i << " <= \"01\";" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "else" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << tab << "readFrame_" << i << " <= \"00\";" << endl
                               << tab << tab << tab << tab << tab << tab << tab << tab << "end if;" << endl
                               << tab << tab << tab << tab << tab << tab << "end case;" << endl
                               << tab << tab << tab << tab << "end case;" << endl
                               << tab << tab << tab << "end if;" << endl
                               << tab << tab << "end if;" << endl
                               << tab << "end if;" << endl
                               << "end process;" << endl << endl;



                }
                else
                {
                    stringstream e;
                    e << "Standard case (2 BRAMs on output) not implemented yet!";
                    THROWERROR(e.str());
                }
            }
        }
    }

    string OutputLayer::getOutputSignalName(int feature)
    {
        return "R_"+to_string(feature);
    }

    string OutputLayer::getInputSignalName(int feature)
    {
        return "X_"+to_string(feature);
    }

    unsigned int OutputLayer::getIndexWordSize() const {
        if(this->myMaxPtr==nullptr) return 0;
        return this->myMaxPtr->getIndexWordSize();
    }

}//namespace flopoco
