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

// include the header of the Operator
#include "ControlledMemory.hpp"
#include "BlockRam.hpp"
#include "MemoryController.hpp"

using namespace std;
namespace flopoco {




    ControlledMemory::ControlledMemory(Target* target, unsigned int dataWidth_, unsigned int addressWidth_, bool outputOverwriteValue, bool needGetNewDataReset, bool hasSpecificAddressOnResetPort, bool hasSpecificAddressOnGetNewDataResetPort) :
        Operator(target), dataWidth(dataWidth_), addressWidth(addressWidth_) {

        this->useNumericStd();

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="ControlledMemory";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // definition of the name of the operator
        ostringstream name;
        name << "ControlledMemory_dataWidth_" << dataWidth << "_addressWidth_" << addressWidth << ((outputOverwriteValue==true)?("_1"):("_0")) << ((needGetNewDataReset==true)?("_1"):("_0")) << ((hasSpecificAddressOnResetPort==true)?("_1"):("_0")) << ((hasSpecificAddressOnGetNewDataResetPort==true)?("_1"):("_0"));
        setName(name.str());

        // add in/out
        addInput("newStep",1);
        addInput("validData_i",1);
        addInput("getNewData",1);
        addInput("data_i",dataWidth);
        if(outputOverwriteValue==true)
        {
            addOutput("data_overwrite_o",dataWidth);
            addInput("getNewData_overwrite",1);
            addInput("getNewData_overwrite_reset",1);
        }

        if(hasSpecificAddressOnResetPort==true)
        {
            addInput("reset_address",this->addressWidth);
        }

        if(needGetNewDataReset==true)
        {
            addInput("getNewData_reset",1);
        }

        if(hasSpecificAddressOnGetNewDataResetPort==true)
        {
            addInput("getNewData_reset_address",this->addressWidth);
        }

        addOutput("validData_o",1);
        addOutput("data_o",dataWidth);

        // signals between BRAM and memory controller
        declare("DataWrite",dataWidth);
        declare("AddressWrite",addressWidth+1);
        declare("WriteEnable",1);
        declare("DataRead",dataWidth);
        declare("AddressRead",addressWidth+1);

        // BRAM
        BlockRam* bramOp = new BlockRam(target, dataWidth, addressWidth+1,outputOverwriteValue);
        this->addSubComponent(bramOp);
        inPortMap(bramOp, "Address_W_in", "AddressWrite");
        inPortMap(bramOp, "Address_R_in", "AddressRead");
        inPortMap(bramOp, "Data_in", "DataWrite");
        inPortMap(bramOp, "WriteEnable", "WriteEnable");
        outPortMap(bramOp, "Data_out", "DataRead", false);
        if(outputOverwriteValue==true)
        {
            outPortMap(bramOp, "Data2_out", "DataRead_2", true);

            declare("Address2Read",addressWidth+1);
            inPortMap(bramOp,"Address_R2_in","Address2Read");
        }
        vhdl << instance(bramOp, "BRAM_instance");

        // memory controller
        MemoryController* memOp = new MemoryController(target, dataWidth, addressWidth,outputOverwriteValue,needGetNewDataReset,hasSpecificAddressOnResetPort,hasSpecificAddressOnGetNewDataResetPort);
        this->addSubComponent(memOp);
        inPortMap(memOp,"newStep","newStep");
        inPortMap(memOp,"validData_i","validData_i");
        inPortMap(memOp,"getNewData","getNewData");
        inPortMap(memOp,"data_i","data_i");
        inPortMap(memOp,"bramDataRead","DataRead");

        outPortMap(memOp,"validData_o","validData_o",false);
        outPortMap(memOp,"data_o","data_o",false);
        outPortMap(memOp,"bramDataWrite","DataWrite",false);
        outPortMap(memOp,"bramAddressRead","AddressRead",false);
        outPortMap(memOp,"bramAddressWrite","AddressWrite",false);
        outPortMap(memOp,"bramWriteEnable","WriteEnable",false);

        if(outputOverwriteValue==true)
        {
            outPortMap(memOp,"data2_o","data_overwrite_o",false);
            inPortMap(memOp,"bramData2Read","DataRead_2");
            outPortMap(memOp,"bramAddress2Read","Address2Read",false);
            inPortMap(memOp,"getNewData2","getNewData_overwrite");
            inPortMap(memOp,"getNewData2_reset","getNewData_overwrite_reset");
        }

        if(needGetNewDataReset==true)
        {
            inPortMap(memOp,"getNewData_reset","getNewData_reset");
        }

        if(hasSpecificAddressOnGetNewDataResetPort==true)
        {
            inPortMap(memOp,"getNewData_reset_address","getNewData_reset_address");
        }

        if(hasSpecificAddressOnResetPort==true)
        {
            inPortMap(memOp,"reset_address","reset_address");
        }

        vhdl << instance(memOp, "MemoryController_instance");
    }

}//namespace flopoco
