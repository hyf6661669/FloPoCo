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




    ControlledMemory::ControlledMemory(Target* target, unsigned int dataWidth_, unsigned int addressWidth_) :
        Operator(target), dataWidth(dataWidth_), addressWidth(addressWidth_) {

        this->useNumericStd();

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="ControlledMemory";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // definition of the name of the operator
        ostringstream name;
        name << "ControlledMemory_dataWidth_" << dataWidth << "_addressWidth_" << addressWidth;
        setName(name.str());

        // add in/out
        addInput("newStep",1);
        addInput("validData_i",1);
        addInput("getNewData",1);
        addInput("data_i",dataWidth);

        addOutput("validData_o",1);
        addOutput("data_o",dataWidth);

        // signals between BRAM and memory controller
        declare("DataWrite",dataWidth);
        declare("AddressWrite",addressWidth+1);
        declare("WriteEnable",1);
        declare("DataRead",dataWidth);
        declare("AddressRead",addressWidth+1);

        // BRAM
        BlockRam* bramOp = new BlockRam(target, dataWidth, addressWidth+1);
        this->addSubComponent(bramOp);
        inPortMap(bramOp, "Address_W_in", "AddressWrite");
        inPortMap(bramOp, "Address_R_in", "AddressRead");
        inPortMap(bramOp, "Data_in", "DataWrite");
        inPortMap(bramOp, "WriteEnable", "WriteEnable");
        outPortMap(bramOp, "Data_out", "DataRead", false);
        vhdl << instance(bramOp, "BRAM_instance");

        // memory controller
        MemoryController* memOp = new MemoryController(target, dataWidth, addressWidth);
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
        vhdl << instance(memOp, "MemoryController_instance");
    }

}//namespace flopoco
