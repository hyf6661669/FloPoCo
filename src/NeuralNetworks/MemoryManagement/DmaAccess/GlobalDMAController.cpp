// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "GlobalDMAController.hpp"

#include "NeuralNetworks/Layers/Layer.hpp"
#include "NeuralNetworks/Utility/DataGuard.hpp"
#include "PrimitiveComponents/GenericLut.hpp"
#include "PrimitiveComponents/GenericMux.hpp"

using namespace std;
namespace flopoco {

    GlobalDMAController::GlobalDMAController(Target* target, unsigned int numberOfInputs_) :
            Operator(target), numberOfInputs(numberOfInputs_) {

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="GlobalDMAController";

        // use numeric_std
        useNumericStd();

        // author
        setCopyrightString("Nicolai Fiege, 2018");

        // definition of the name of the operator
        ostringstream name;
        name << "GlobalDMAController_" << this->numberOfInputs;
        setName(name.str());

        /////////////////////////////////////
        // communication with other layers //
        /////////////////////////////////////
        for(unsigned int i = 0; i < this->numberOfInputs; i++)
        {
            // actual data
            addOutput(this->getDataPortName(i),32);
			// last data
			addOutput(this->getLastDataPortName(i),1);
            // this layer gets valid data
            addOutput(Layer::getValidDataName(i),1);
            // this layer can accept new data in the next clock cycle
            //addInput(Layer::getGetNewDataName(i),1);
            // this layer wants to start a new read-access
            addInput(this->getNewReadAccessPortName(i),1);
            // start address this layer wants to read weights from
            addInput(this->getNewStartAddressPortName(i),32);
            // number of bytes this layer wants to read
            addInput(this->getNumberOfTotalBytesPortName(i),23);

            // std_logic-signal for Priority-LUT
            this->vhdl << this->declare(this->getNewReadAccessPortName(i)+"_std",1,false) << " <= " << this->getNewReadAccessPortName(i) << "(0);" << endl;
        }

        ////////////////////////////
        // communication with DMA //
        ////////////////////////////

        // ready to get new weight from DMA
        addOutput(GlobalDMAController::getDMAReadyPortName(),1);
        // DMA sends valid data
        addInput(GlobalDMAController::getDMAValidPortName(),1);
        // DMA sends last data from this packet
        addInput(GlobalDMAController::getDMALastPortName(),1);
        // actual data
        addInput(GlobalDMAController::getDMADataPortName(),32);

        // starts new DMA transfer on rising edge
        addOutput(GlobalDMAController::getNewDMAAccessPortName(),1);
        // start address to read data from
        addOutput(GlobalDMAController::getDMAStartAddressPortName(),32);
        // number of bytes to read
        addOutput(GlobalDMAController::getDMANumberOfBytesPortName(),23);
        // DMA is programmed for a new read access
        addInput(GlobalDMAController::getDMAFinishedProgrammingPortName(),1);
        // DMA is ready to read another chunk of data
        addInput(GlobalDMAController::getDMAReadyToBeProgrammedPortName(),1);

        ////////////////////////////////
        // keep track of new requests //
        ////////////////////////////////
        this->declare("New_Request_On_One_Input",1);
        this->vhdl << "New_Request_On_One_Input <= " << this->getNewReadAccessPortName(0);
        for(unsigned int i=1; i<this->numberOfInputs; i++)
        {
            this->vhdl << " or " << this->getNewReadAccessPortName(i);
        }
        this->vhdl << ";" << endl;

        ////////////////////////////////////////////////////////////
        // MUXs and LUT to select the input with highest priority //
        ////////////////////////////////////////////////////////////
        if(this->numberOfInputs > 1)
        {
            this->declare("Priority_Port_No",ceil(log2(this->numberOfInputs)));
            this->declare("Priority_Port_No_reg",ceil(log2(this->numberOfInputs)));
            // LUT to generate mux-select
            map<unsigned int, unsigned int> LUTMap;
            LUTMap[0] = 0;
            for(unsigned int i=1; i<pow(2,this->numberOfInputs); i++)
            {
                unsigned int reversed = ceil(log2(i+1));
                LUTMap[i] = this->numberOfInputs-reversed;
            }
            GenericLut* LUT = new GenericLut(target,"Priority_Lut",LUTMap,this->numberOfInputs,ceil(log2(this->numberOfInputs)));
            this->addSubComponent(LUT);
            // connect inputs in reversed order
            for(unsigned int i = 0; i < this->numberOfInputs; i++)
            {
                this->inPortMap(LUT,join("i",i),this->getNewReadAccessPortName(this->numberOfInputs-i-1)+"_std");
            }
            for(unsigned int i = 0; i < ceil(log2(this->numberOfInputs)); i++)
            {
                this->outPortMap(LUT,join("o",i),"Priority_Port_No_"+to_string(i), true);
                this->vhdl << "Priority_Port_No(" << i << ") <= Priority_Port_No_" << i << ";" << endl;
            }
            this->vhdl << instance(LUT, "Priority_Lut_instance");

            this->vhdl << "process(clk)" << endl;
            this->vhdl << "begin" << endl;
            this->vhdl << tab << "if(rising_edge(clk)) then" << endl;
            this->vhdl << tab << tab << "if(rst='1') then" << endl;
            this->vhdl << tab << tab << tab << "Priority_Port_No_reg <= (others => '0');" << endl;
            this->vhdl << tab << tab << "elsif(Controller_State=\"00\" and New_Request_On_One_Input = \"1\") then" << endl;
            this->vhdl << tab << tab << tab << "Priority_Port_No_reg <= Priority_Port_No;" << endl;
            this->vhdl << tab << tab << "end if;" << endl;
            this->vhdl << tab << "end if;" << endl;
            this->vhdl << "end process;" << endl;
        }

        if(this->numberOfInputs > 1)
        {
            // next start address mux
            GenericMux* nextStartAddressMux = new GenericMux(target,32,this->numberOfInputs);
            this->addSubComponent(nextStartAddressMux);
            this->inPortMap(nextStartAddressMux,nextStartAddressMux->getSelectName(),"Priority_Port_No_reg");
            for(unsigned int i=0; i<this->numberOfInputs; i++)
            {
                this->inPortMap(nextStartAddressMux,nextStartAddressMux->getInputName(i),this->getNewStartAddressPortName(i));
            }
            this->outPortMap(nextStartAddressMux,nextStartAddressMux->getOutputName(),"Next_Start_Address",true);
            this->vhdl << instance(nextStartAddressMux, "Start_Address_Mux_instance");

            // next number of bytes mux
            GenericMux* nextNumberOfBytesMux = new GenericMux(target,23,this->numberOfInputs);
            this->addSubComponent(nextNumberOfBytesMux);
            this->inPortMap(nextNumberOfBytesMux,nextNumberOfBytesMux->getSelectName(),"Priority_Port_No_reg");
            for(unsigned int i=0; i<this->numberOfInputs; i++)
            {
                this->inPortMap(nextNumberOfBytesMux,nextNumberOfBytesMux->getInputName(i),this->getNumberOfTotalBytesPortName(i));
            }
            this->outPortMap(nextNumberOfBytesMux,nextNumberOfBytesMux->getOutputName(),"Next_Number_Of_Bytes",true);
            this->vhdl << instance(nextNumberOfBytesMux, "Number_Of_Bytes_Mux_instance");
        }
        else
        {
            this->vhdl << declare("Next_Start_Address",32) << " <= " << this->getNewStartAddressPortName(0) << ";" << endl;
            this->vhdl << declare("Next_Number_Of_Bytes",23) << " <= " << this->getNumberOfTotalBytesPortName(0) << ";" << endl;
        }

        ////////////////////////////////////////////////
        // Create Guards for valid data stream signal //
        ////////////////////////////////////////////////
        if(this->numberOfInputs > 1)
        {
            for(unsigned int i = 0; i < this->numberOfInputs; i++)
            {
                DataGuard* dat = new DataGuard(target,1,ceil(log2(this->numberOfInputs)),i);
                addSubComponent(dat);
                inPortMap(dat,"Data_in",GlobalDMAController::getDMAValidPortName());
                inPortMap(dat,"Guard_in","Priority_Port_No_reg");
                outPortMap(dat,"Data_out",Layer::getValidDataName(i),false);
                this->vhdl << instance(dat, "ValidDataGuard_"+to_string(i));
            }
        }
        else
        {
            this->vhdl << Layer::getValidDataName(0) << " <= " << this->getDMAValidPortName() << ";" << endl;
        }

        ///////////////////////////////////////////////
        // Create Guards for last data stream signal //
        ///////////////////////////////////////////////
        if(this->numberOfInputs > 1)
        {
            for(unsigned int i = 0; i < this->numberOfInputs; i++)
            {
                DataGuard* dat = new DataGuard(target,1,ceil(log2(this->numberOfInputs)),i);
                addSubComponent(dat);
                inPortMap(dat,"Data_in",GlobalDMAController::getDMALastPortName());
                inPortMap(dat,"Guard_in","Priority_Port_No_reg");
                outPortMap(dat,"Data_out",this->getLastDataPortName(i),false);
                this->vhdl << instance(dat, "LastDataGuard_"+to_string(i));
            }
        }
        else
        {
            this->vhdl << this->getLastDataPortName(0) << " <= " << GlobalDMAController::getDMALastPortName() << ";" << endl;
        }


        /////////////////////////////////////////////////////////////////
        // FSM to manage wait for new access/program DMA/read from DMA //
        /////////////////////////////////////////////////////////////////
        declare("Controller_State",2);
        // 00: wait for new access
        // 01: init DMA programming and wait until it's programmed
        // 10: read data from DMA
        // 11: wait until DMA can be programmed again
        this->vhdl << "process(clk)" << endl;
        this->vhdl << "begin" << endl;
        this->vhdl << tab << "if(rising_edge(clk)) then" << endl;
        this->vhdl << tab << tab << "if(rst='1') then" << endl;
        this->vhdl << tab << tab << tab << "Controller_State <= \"11\";" << endl;
        this->vhdl << tab << tab << "else" << endl;
        this->vhdl << tab << tab << tab << "case Controller_State is" << endl;
        this->vhdl << tab << tab << tab << tab << "when \"00\" => " << endl;
        this->vhdl << tab << tab << tab << tab << tab << "if(New_Request_On_One_Input = \"1\") then" << endl;
        this->vhdl << tab << tab << tab << tab << tab << tab << "Controller_State <= \"01\";" << endl;
        this->vhdl << tab << tab << tab << tab << tab << "else" << endl;
        this->vhdl << tab << tab << tab << tab << tab << tab << "Controller_State <= \"00\";" << endl;
        this->vhdl << tab << tab << tab << tab << tab << "end if;" << endl;
        this->vhdl << tab << tab << tab << tab << "when \"01\" => " << endl;
        this->vhdl << tab << tab << tab << tab << tab << "if(" << GlobalDMAController::getDMAFinishedProgrammingPortName() << " = \"1\") then" << endl;
        this->vhdl << tab << tab << tab << tab << tab << tab << "Controller_State <= \"10\";" << endl;
        this->vhdl << tab << tab << tab << tab << tab << "else" << endl;
        this->vhdl << tab << tab << tab << tab << tab << tab << "Controller_State <= \"01\";" << endl;
        this->vhdl << tab << tab << tab << tab << tab << "end if;" << endl;
        this->vhdl << tab << tab << tab << tab << "when \"10\" => " << endl;
        this->vhdl << tab << tab << tab << tab << tab << "if(" << GlobalDMAController::getDMALastPortName() << " = \"1\") then" << endl;
        this->vhdl << tab << tab << tab << tab << tab << tab << "Controller_State <= \"11\";" << endl;
        this->vhdl << tab << tab << tab << tab << tab << "else" << endl;
        this->vhdl << tab << tab << tab << tab << tab << tab << "Controller_State <= \"10\";" << endl;
        this->vhdl << tab << tab << tab << tab << tab << "end if;" << endl;
        this->vhdl << tab << tab << tab << tab << "when \"11\" => " << endl;
        this->vhdl << tab << tab << tab << tab << tab << "if(" << GlobalDMAController::getDMAReadyToBeProgrammedPortName() << " = \"1\") then" << endl;
        this->vhdl << tab << tab << tab << tab << tab << tab << "Controller_State <= \"00\";" << endl;
        this->vhdl << tab << tab << tab << tab << tab << "else" << endl;
        this->vhdl << tab << tab << tab << tab << tab << tab << "Controller_State <= \"11\";" << endl;
        this->vhdl << tab << tab << tab << tab << tab << "end if;" << endl;
        this->vhdl << tab << tab << tab << tab << "when others => " << endl;
        this->vhdl << tab << tab << tab << tab << tab << "Controller_State <= \"11\";" << endl;
        this->vhdl << tab << tab << tab << "end case;" << endl;
        this->vhdl << tab << tab << "end if;" << endl;
        this->vhdl << tab << "end if;" << endl;
        this->vhdl << "end process;" << endl;

        // handle dma programming
        this->vhdl << "process(Controller_State)" << endl;
        this->vhdl << "begin" << endl;
        this->vhdl << tab << "case Controller_State is" << endl;
        this->vhdl << tab << tab << "when \"01\" => " << endl;
        this->vhdl << tab << tab << tab << GlobalDMAController::getNewDMAAccessPortName() << " <= \"1\";" << endl;
        this->vhdl << tab << tab << "when others => " << endl;
        this->vhdl << tab << tab << tab << GlobalDMAController::getNewDMAAccessPortName() << " <= \"0\";" << endl;
        this->vhdl << tab << "end case;" << endl;
        this->vhdl << "end process;" << endl;

        this->vhdl << GlobalDMAController::getDMANumberOfBytesPortName() << " <= Next_Number_Of_Bytes;" << endl;
        this->vhdl << GlobalDMAController::getDMAStartAddressPortName() << " <= Next_Start_Address;" << endl;

        // handle reading from DMA
        this->vhdl << "process(Controller_State)" << endl;
        this->vhdl << "begin" << endl;
        this->vhdl << tab << "case Controller_State is" << endl;
        this->vhdl << tab << tab << "when \"10\" => " << endl;
        this->vhdl << tab << tab << tab << GlobalDMAController::getDMAReadyPortName() << " <= \"1\";" << endl;
        this->vhdl << tab << tab << "when others => " << endl;
        this->vhdl << tab << tab << tab << GlobalDMAController::getDMAReadyPortName() << " <= \"0\";" << endl;
        this->vhdl << tab << "end case;" << endl;
        this->vhdl << "end process;" << endl;

        ///////////////////////////////////////////////////////////////
        // Connect all Data-to-Layer-signals with data stream signal //
        ///////////////////////////////////////////////////////////////
        for(unsigned int i=0; i<this->numberOfInputs; i++)
        {
            this->vhdl << this->getDataPortName(i) << " <= " << GlobalDMAController::getDMADataPortName() << ";" << endl;
        }

    }

    string GlobalDMAController::getNewReadAccessPortName(unsigned int portNumber_) const
    {
        if(portNumber_ > this->numberOfInputs)
        {
            THROWERROR("Requested port doesn't exits");
        }
        return "newReadAccess_"+to_string(portNumber_);
    }

    string GlobalDMAController::getDataPortName(unsigned int portNumber_) const
    {
        if(portNumber_ > this->numberOfInputs)
        {
            THROWERROR("Requested port doesn't exits");
        }
        return "Data_"+to_string(portNumber_);
    }

    string GlobalDMAController::getNewStartAddressPortName(unsigned int portNumber_) const
    {
        if(portNumber_ > this->numberOfInputs)
        {
            THROWERROR("Requested port doesn't exits");
        }
        return "newStartAddress_"+to_string(portNumber_);
    }

    string GlobalDMAController::getLastDataPortName(unsigned int portNumber_) const
    {
        if(portNumber_ > this->numberOfInputs)
        {
            THROWERROR("Requested port doesn't exits");
        }
        return "lastData_"+to_string(portNumber_);
    }

    string GlobalDMAController::getNumberOfTotalBytesPortName(unsigned int portNumber_) const
    {
        if(portNumber_ > this->numberOfInputs)
        {
            THROWERROR("Requested port doesn't exits");
        }
        return "numberOfBytes_"+to_string(portNumber_);
    }

}//namespace flopoco